#include <iostream>
#include <cuda_runtime.h>

__global__ void helloFromGPU()
{
    printf("Hello from GPU thread %d\n", threadIdx.x);
}

int main()
{
    helloFromGPU<<<1, 4>>>();

    cudaError_t error = cudaDeviceSynchronize();

    if (error != cudaSuccess)
    {
        std::cerr << "CUDA error: "
                  << cudaGetErrorString(error)
                  << '\n';

        return 1;
    }

    std::cout << "CUDA program completed.\n";

    return 0;
}