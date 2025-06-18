#ifndef MODEL_WRAPPER_HH
#define MODEL_WRAPPER_HH

#ifdef __cplusplus
extern "C" {
#endif

int predict(float *input_array, int num_elements, int input_size);

#ifdef __cplusplus
}
#endif

#endif // MODEL_WRAPPER_HH
