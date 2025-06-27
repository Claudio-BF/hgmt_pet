#include <filesystem>
#include <omp.h>
#include <torch/script.h>
#include <torch/torch.h>
#include <unistd.h>
using namespace std;
static string model_path() {
  char buf[1024];
  ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  buf[len] = '\0';
  string exePath(buf);
  string dirPath = exePath.substr(0, exePath.find_last_of('/') + 1);
  string modelPath = dirPath + "neural_network/chooser.pt";
  return modelPath;
}
// Configure threading before any PyTorch operations
static void configure_single_threading() {
  // Set environment variables first (most reliable)
  setenv("OMP_NUM_THREADS", "1", 1);
  setenv("MKL_NUM_THREADS", "1", 1);

  // Then set PyTorch threading
  at::set_num_interop_threads(1);
  at::set_num_threads(1);
  torch::set_num_threads(1);
  omp_set_num_threads(1);
  omp_set_dynamic(0);
}

// Internal C++ implementation
static int predict_impl(float *input_array, int num_elements, int input_size) {
  thread_local torch::jit::script::Module module = []() {
    auto m = torch::jit::load(model_path());
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
  // Call threading configuration once at the very beginning
  static once_flag threading_configured;
  call_once(threading_configured, configure_single_threading);

  return predict_impl(input_array, num_elements, input_size);
}
}
