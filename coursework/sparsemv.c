#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <immintrin.h>  // 包含 AVX, AVX2 指令
#include "sparsemv.h"

/**
 * @brief Compute matrix vector product (y = A*x)
 * 
 * @param A Known matrix
 * @param x Known vector
 * @param y Return vector
 * @return int 0 if no error
 */
int sparsemv(struct mesh *A, const double * const x, double * const y)
{

  const int nrow = (const int) A->local_nrow;

  // #pragma omp parallel for
  for (int i=0; i< nrow; i++) {
      double sum = 0.0;
      const double * const cur_vals = (const double * const) A->ptr_to_vals_in_row[i];
      const int * const cur_inds = (const int * const) A->ptr_to_inds_in_row[i];
      const int cur_nnz = (const int) A->nnz_in_row[i];

      
      for (int j=0; j< cur_nnz; j++) {
        sum += cur_vals[j]*x[cur_inds[j]];
      }
      y[i] = sum;

      // __m256d vec_sum = _mm256_setzero_pd();  // 初始化 4 个 double 累加和
      // int j = 0;

      // // **SIMD 处理 4 个元素**
      // for (; j + 3 < cur_nnz; j += 4) {
      //     __m256d vec_vals = _mm256_loadu_pd(&cur_vals[j]);  // 读取 4 个矩阵值

      //     // 使用 gather 从 x 中加载 4 个不连续的元素
      //     __m128i indices = _mm_loadu_si128((__m128i*)&cur_inds[j]);  
      //     __m256d vec_x = _mm256_i32gather_pd(x, indices, sizeof(double));

      //     vec_sum = _mm256_add_pd(vec_sum, _mm256_mul_pd(vec_vals, vec_x));

      // }

      // // **手动累加 SIMD 结果**
      // double temp[4];
      // _mm256_storeu_pd(temp, vec_sum);
      // sum += temp[0] + temp[1] + temp[2] + temp[3];

      // // **处理剩余的非 4 的倍数**
      // for (; j < cur_nnz; j++) {
      //     sum += cur_vals[j] * x[cur_inds[j]];
      // }

      // y[i] = sum;





    }
  return 0;
}
