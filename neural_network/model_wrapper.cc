#include <omp.h>
#include <torch/script.h>
#include <torch/torch.h>

// Add this before your predict_impl function
static void configure_single_threading() {
  at::set_num_interop_threads(1);
  at::set_num_threads(1);
  torch::set_num_threads(1);
  omp_set_num_threads(1);
  omp_set_dynamic(0);
}

const std::string MODEL_PATH =
    "/home/claud/Work/hgmt_pet_scanner/neural_network/chooser.pt";

// Internal C++ implementation
static int predict_impl(float *input_array, int num_elements, int input_size) {
  thread_local torch::jit::script::Module module = []() {
    configure_single_threading();
    auto m = torch::jit::load(MODEL_PATH);
    m.eval();
    return m;
  }();

  torch::Tensor input =
      torch::from_blob(input_array, {1, num_elements, input_size},
                       torch::kFloat)
          .clone();

  torch::Tensor output;
  {
    torch::InferenceMode guard;
    output = module.forward({input}).toTensor();
  }

  return torch::argmax(output).item<int>();
}

// C-compatible wrapper
extern "C" {
int predict(float *input_array, int num_elements, int input_size) {
  return predict_impl(input_array, num_elements, input_size);
}
}
