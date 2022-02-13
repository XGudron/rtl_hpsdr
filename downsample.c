/*	This file downsample.c is part of rtl_hpsdr.
 *
 *	rtl_hpsdr - an RTL to HPSDR software translation server
 *	Copyright (C) 2014 Richard Koch
 *
 *	rtl_hpsdr is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	rtl_hpsdr is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with rtl_hpsdr.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "rtl_hpsdr.h"
#include "biquad.h"

#ifdef INCLUDE_NEON
#define vector_type float32x4_t
#define vector_zero vmovq_n_f32(0)
#define vector_load(x) vld1q_f32(x)
#define vector_mac(x,y,z) vmlaq_f32((x), vld1q_f32((y)), (z))
#define vector_store(x,y,z) vst1q_f32((x), vaddq_f32((y), vcombine_f32(vget_low_f32((z)), vget_high_f32((y)))))
#elif defined INCLUDE_SSE2
#define vector_type __m128
#define vector_zero _mm_setzero_ps()
#define vector_load(x) _mm_load_ps(x)
#define vector_mac(x,y,z) _mm_add_ps((x), _mm_mul_ps(_mm_load_ps((y)), (z)))
#define vector_store(x,y,z) _mm_store_ps((x), _mm_add_ps((y), _mm_shuffle_ps((y), (y), _MM_SHUFFLE(1,0,3,2))))
#else
#define vector_type float
#define vector_load(x) (*x)
#endif

iir_filter_t *filter;
static float32_t IIR_LPF_Coeffs[BIQUAD_COEFF_IN_STAGE * IIR_LPF_STAGES] = {0};

void arm_biquad_cascade_df2T_f32_rolled(const arm_biquad_cascade_df2T_instance_f32 * S,const float32_t * pSrc,float32_t * pDst,uint32_t blockSize)
{
  const float32_t *pIn = pSrc;                         /* Source pointer */
        float32_t *pOut = pDst;                        /* Destination pointer */
        float32_t *pState = S->pState;                 /* State pointer */
  const float32_t *pCoeffs = S->pCoeffs;               /* Coefficient pointer */
        float32_t acc1;                                /* Accumulator */
        float32_t b0, b1, b2, a1, a2;                  /* Filter coefficients */
        float32_t Xn1;                                 /* Temporary input */
        float32_t d1, d2;                              /* State variables */
        uint32_t sample, stage = S->numStages;         /* Loop counters */

  do
  {
     /* Reading the coefficients */
     b0 = pCoeffs[0];
     b1 = pCoeffs[1];
     b2 = pCoeffs[2];
     a1 = pCoeffs[3];
     a2 = pCoeffs[4];

     /* Reading the state values */
     d1 = pState[0];
     d2 = pState[1];

     pCoeffs += 5U;

      /* Initialize blkCnt with number of samples */
      sample = blockSize;

      while (sample > 0U) {
        Xn1 = *pIn++;

        acc1 = b0 * Xn1 + d1;

        d1 = b1 * Xn1 + d2;
        d1 += a1 * acc1;

        d2 = b2 * Xn1;
        d2 += a2 * acc1;

        *pOut++ = acc1;

        /* decrement loop counter */
        sample--;
      }

      /* Store the updated state variables back into the state array */
      pState[0] = d1;
      pState[1] = d2;

      pState += 2U;

      /* The current stage output is given as the input to the next stage */
      pIn = pDst;

      /* Reset the output working pointer */
      pOut = pDst;

      /* decrement loop counter */
      stage--;

   } while (stage > 0U);
}

void fill_biquad_coeffs(iir_filter_t *filter, float32_t *coeffs, uint8_t sect_num)
{
	//transpose and save coefficients
	uint16_t ind = 0;
	for(uint8_t sect = 0; sect < sect_num; sect++)
	{
		coeffs[ind + 0] = filter->b[sect * 3 + 0];
		coeffs[ind + 1] = filter->b[sect * 3 + 1];
		coeffs[ind + 2] = filter->b[sect * 3 + 2];
		coeffs[ind + 3] = -filter->a[sect * 3 + 1];
		coeffs[ind + 4] = -filter->a[sect * 3 + 2];
		ind += 5;
	}
}

void arm_biquad_cascade_df2T_init_f32(
  arm_biquad_cascade_df2T_instance_f32 * S,
  uint8_t numStages,
  float32_t * pCoeffs,
  float32_t * pState)
{
  /* Assign filter stages */
  S->numStages = numStages;

  /* Assign coefficient pointer */
  S->pCoeffs = pCoeffs;

  /* Clear state buffer and size is always 2 * numStages */
  memset(pState, 0, (2u * (uint32_t) numStages) * sizeof(float32_t));

  /* Assign state pointer */
  S->pState = pState;
}

void downsample_init(struct main_cb* mcb) {
	filter = biquad_create(IIR_LPF_STAGES);
	biquad_init_lowpass(filter, RTL_SAMPLE_RATE, mcb->output_rate / 2);
	printf("Init LPF %d\n", mcb->output_rate / 2);
	fill_biquad_coeffs(filter, IIR_LPF_Coeffs, IIR_LPF_STAGES);

	for(uint8_t i = 0; i < mcb->active_num_rcvrs; i++) {
		arm_biquad_cascade_df2T_init_f32(&(mcb->rcb[i].IIR_LPF_I), IIR_LPF_STAGES, IIR_LPF_Coeffs, (float32_t *)&(mcb->rcb[i].IIR_LPF_I_State[0]));
		arm_biquad_cascade_df2T_init_f32(&(mcb->rcb[i].IIR_LPF_Q), IIR_LPF_STAGES, IIR_LPF_Coeffs, (float32_t *)&(mcb->rcb[i].IIR_LPF_Q_State[0]));
		mcb->rcb[i].LPF_inited = true;
	}
}

void decimate(float32_t *buff, uint32_t size, uint8_t rate) {
	for(uint32_t i = 0; i < size / rate; i++){
		buff[i] = buff[i * rate];
	}
}

void
downsample(struct rcvr_cb * rcb) {
	if(!rcb->LPF_inited)
		downsample_init(rcb->mcb);

	arm_biquad_cascade_df2T_f32_rolled(&(rcb->IIR_LPF_I), &(rcb->iq_buf_I[0]), &(rcb->iq_buf_I[0]), RTL_READ_COUNT / 2);
	arm_biquad_cascade_df2T_f32_rolled(&(rcb->IIR_LPF_Q), &(rcb->iq_buf_Q[0]), &(rcb->iq_buf_Q[0]), RTL_READ_COUNT / 2);

	uint8_t decim_rate = DOWNSAMPLE_192;
	if(rcb->mcb->output_rate == 48000)
		decim_rate = DOWNSAMPLE_192 * 4;
	if(rcb->mcb->output_rate == 96000)
		decim_rate = DOWNSAMPLE_192 * 2;
	if(rcb->mcb->output_rate == 384000)
		decim_rate = DOWNSAMPLE_192 / 2;

	decimate(&(rcb->iq_buf_I[0]), RTL_READ_COUNT / 2, decim_rate);
	decimate(&(rcb->iq_buf_Q[0]), RTL_READ_COUNT / 2, decim_rate);

	uint32_t j = 0;
	float *out = &(rcb->iqSamples[rcb->iqSamples_remaining * 2]);
	for(uint32_t i = 0; i < (RTL_READ_COUNT / 2 / decim_rate); i++)
	{
		static const float A1 = (1.0f - powf(2, -7)); // (1-2^(-11))

		//I correct + swap + copy
		static float I_x_prev = 0.0f;
		static float I_y_prev = 0.0f;
		float sampleIn = rcb->iq_buf_I[i];
		float sampleOut = 0;
		float delta_x = sampleIn - I_x_prev;
		float a1_y_prev = A1 * I_y_prev;
		sampleOut = delta_x + a1_y_prev;
		I_x_prev = sampleIn;
		I_y_prev = sampleOut;
		out[j+1] = sampleOut;

		//Q correct + swap + copy
		static float Q_x_prev = 0.0f;
		static float Q_y_prev = 0.0f;
		sampleIn = rcb->iq_buf_Q[i];
		sampleOut = 0;
		delta_x = sampleIn - Q_x_prev;
		a1_y_prev = A1 * Q_y_prev;
		sampleOut = delta_x + a1_y_prev;
		Q_x_prev = sampleIn;
		Q_y_prev = sampleOut;
		out[j] = sampleOut;

		//iq step
		j=j+2;
	}
	rcb->LPF_downsampled = true;
}
