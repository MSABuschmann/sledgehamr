#include "pytorch_integrator.h"

namespace sledgehamr {

IntegratorPytorch::IntegratorPytorch(Sledgehamr *owner) : Integrator{owner} {
    ParseParams();
    LoadPytorchModel();
}

void IntegratorPytorch::ParseParams() {
    amrex::ParmParse pp("integrator");
    pp.get("pytorch_model_filename", pytorch_model_filename);
}

void IntegratorPytorch::LoadPytorchModel() {
    try {
        module = torch::jit::load(pytorch_model_filename);
    } catch (const c10::Error &e) {
        amrex::Abort("Error loading the Pytorch model!\n");
    }

#ifdef AMREX_USE_CUDA
    torch::Device device0(torch::kCUDA);
    module.to(device0);
    amrex::Print() << "Copying model to GPU." << std::endl;

    // set tensor options
    tensoropt = torch::TensorOptions().dtype(dtype0).device(device0);
#else
    tensoropt = torch::TensorOptions().dtype(dtype0);
#endif
}

void IntegratorPytorch::Integrate(LevelData &mf_old, LevelData &mf_new,
                                  const int lev, const double dt,
                                  const double dx) {
    const int ncomp_in = 1;
    const int ncomp_out = 2;

    for (amrex::MFIter mfi(mf_new); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.validbox();

        const amrex::Array4<double> &in_arr = mf_new.array(mfi);
        const amrex::Array4<double> &out_arr = mf_old.array(mfi);

        const amrex::IntVect bx_lo = bx.smallEnd();
        const amrex::IntVect nbox = bx.size();
        int ncell = nbox[0] * nbox[1] * nbox[2];

        amrex::Gpu::ManagedVector<double> aux(ncell * ncomp_in);
        double *AMREX_RESTRICT auxPtr = aux.dataPtr();

        amrex::ParallelFor(
            bx, ncomp_in,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                int ii = i - bx_lo[0];
                int jj = j - bx_lo[1];
                int kk = k - bx_lo[2];
                int index = kk * nbox[0] * nbox[1] + jj * nbox[0] + ii;

                // array order is row-based [index][comp]
                auxPtr[index * ncomp_in + n] = in_arr(i, j, k, n);
            });

        at::Tensor inputs_torch =
            torch::from_blob(auxPtr, {ncell, ncomp_in}, tensoropt);

        at::Tensor outputs_torch = module.forward({inputs_torch}).toTensor();
        outputs_torch = outputs_torch.to(dtype0);

#ifdef AMREX_USE_CUDA
        auto outputs_torch_acc = outputs_torch.packed_accessor64<double, 2>();
#else
        auto outputs_torch_acc = outputs_torch.accessor<double, 2>();
#endif

        amrex::ParallelFor(
            bx, ncomp_out,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                int ii = i - bx_lo[0];
                int jj = j - bx_lo[1];
                int kk = k - bx_lo[2];
                int index = kk * nbox[0] * nbox[1] + jj * nbox[0] + ii;
                out_arr(i, j, k, n) = outputs_torch_acc[index][n];
            });
    }
}

}; // namespace sledgehamr
