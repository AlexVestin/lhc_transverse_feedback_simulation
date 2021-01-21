import PyHEADTAIL_feedback
from PyHEADTAIL_feedback.feedback import OneboxFeedback
from PyHEADTAIL_feedback.processors.multiplication import ChargeWeighter
from PyHEADTAIL_feedback.processors.register import TurnFIRFilter
from PyHEADTAIL_feedback.processors.convolution import Lowpass, FIRFilter
from PyHEADTAIL_feedback.processors.resampling import DAC, HarmonicADC, BackToOriginalBins, Upsampler
from PyHEADTAIL_feedback.MD4063_filter_functions import calculate_coefficients_3_tap, calculate_hilbert_notch_coefficients

from PyHEADTAIL_feedback.processors.misc import Bypass
from PyHEADTAIL.particles.slicing import UniformBinSlicer

from scipy.constants import c

import numpy as np

class Damper:
	def __init__(self,machine):
		self.dampingtime = 20.
		self.gain = 2./self.dampingtime

		#lowpass100kHz = [1703, 1169, 1550, 1998, 2517, 3108, 3773, 4513, 5328, 6217, 7174, 8198, 9282, 10417, 11598, 12813, 14052, 15304, 16555, 17793, 19005, 20176, 21294, 22345, 23315, 24193, 24969, 25631, 26171, 26583, 26860, 27000, 27000, 26860, 26583, 26171, 25631, 24969, 24193, 23315, 22345, 21294, 20176, 19005, 17793, 16555, 15304, 14052, 12813, 11598, 10417, 9282, 8198, 7174, 6217, 5328, 4513, 3773, 3108, 2517, 1998, 1550, 1169, 1703]

		#lowpassEnhanced = [490,177,-478,-820,-370,573,1065,428,-909, -1632,-799,1015, 2015,901,-1592,-3053,-1675,1642, 3670,1841,-2828,-6010,-3929,2459,7233,4322,-6384,-17305,-18296,-5077,16097,32000, 32000,16097,-5077,-18296,-17305,-6384,4322, 7233,2459,-3929,-6010,-2828,1841,3670,1642,-1675,-3053,-1592,901,2015,1015, -799,-1632,-909,428,1065,573,-370,-820,-478,177,490]

		self.lowpass20MHz = [38,118,182,112,-133,-389,-385,-45,318,257,-259,-665,-361,473,877,180,-996,-1187,162,1670,1329,-954, -2648, -1219,2427,4007,419,-5623, -6590,2893,19575,32700,32700,19575, 2893,-6590,-5623,419,4007,2427,-1219,-2648, -954, 1329,1670, 162,-1187,-996,180,877,473,-361,-665,-259, 257,318,-45,-385,-389,-133,
		112,182,118,38]

		self.phaseEqualizer = [2,4,7,10,12,16,19,22,27,31,36,42,49,57,67,77,90,104,121,141,164,191,223,261,305, 358,422,498,589,700,836,1004,1215,1483,1832,2301, 2956,3944,5600,9184,25000,-16746,-4256,-2056,-1195,-769,-523,-372,-271,-202,-153, -118,-91,-71,-56,-44,-34,-27,-20,-15,-11,-7,-4,-1] 
		    
		self.FIR_phase_filter = np.array(self.phaseEqualizer)    
		self.FIR_phase_filter = self.FIR_phase_filter/float(np.sum(self.FIR_phase_filter))

		self.FIR_gain_filter = np.array(self.lowpass20MHz)
		self.FIR_gain_filter = self.FIR_gain_filter/float(np.sum(self.lowpass20MHz))


		 # Cut-off frequency of the kicker system
		self.fc=1.0e6
		self.ADC_bits = 16 
		self.ADC_range = (-1e-3, 1e-3)

		# signal processing delay in turns before the first measurements is applied
		self.delay = 1
		self.extra_adc_bins = 10
		# betatron phase advance between the pickup and the kicker. The value 0.25 
		# corresponds to the 90 deg phase change from from the pickup measurements
		# in x-plane to correction kicks in xp-plane.

		self.additional_phase = 0.25 # Kicker-to-pickup phase advance 0 deg
		#    additional_phase = 0. # Kicker-to-pickup phase advance 90 deg


		self.f_RF = 1./(machine.circumference/c/(float(machine.h_RF)))
		self.turn_phase_filter_x = calculate_hilbert_notch_coefficients(machine.Q_x, self.delay, self.additional_phase)
		self.turn_phase_filter_y = calculate_hilbert_notch_coefficients(machine.Q_y, self.delay, self.additional_phase)

		#turn_phase_filter_x = calculate_coefficients_3_tap(machine.Q_x, delay, additional_phase)
		#turn_phase_filter_y = calculate_coefficients_3_tap(machine.Q_y, delay, additional_phase)

		print('f_RF: ' + str(self.f_RF))


		self.processors_detailed_x = [
		        Bypass(),
		        ChargeWeighter(normalization = 'segment_average'),
		#         NoiseGenerator(RMS_noise_level, debug=False),
		        HarmonicADC(1*self.f_RF/10., self.ADC_bits, self.ADC_range,
		                    n_extras=self.extra_adc_bins),
		        TurnFIRFilter(self.turn_phase_filter_x, machine.Q_x, delay=self.delay),
		        FIRFilter(self.FIR_phase_filter, zero_tap = 40),
		        Upsampler(3, [1.5,1.5,0]),
		        FIRFilter(self.FIR_gain_filter, zero_tap = 34),
		        DAC(self.ADC_bits, self.ADC_range),
		        Lowpass(self.fc, f_cutoff_2nd=10*self.fc),
		        BackToOriginalBins(),
		]

		self.processors_detailed_y = [
		        Bypass(),
		        ChargeWeighter(normalization = 'segment_average'),
		#         NoiseGenerator(RMS_noise_level, debug=False),
		        HarmonicADC(1*self.f_RF/10., self.ADC_bits, self.ADC_range,
		                    n_extras=self.extra_adc_bins),
		        TurnFIRFilter(self.turn_phase_filter_y, machine.Q_y, delay = self.delay),
		        FIRFilter(self.FIR_phase_filter, zero_tap = 40),
		        Upsampler(3, [1.5,1.5,0]),
		        FIRFilter(self.FIR_gain_filter, zero_tap = 34),
		        DAC(self.ADC_bits, self.ADC_range),
		        Lowpass(self.fc, f_cutoff_2nd=10*self.fc),
		        BackToOriginalBins(),
		]

		self.sigma_z = 1.2e-9 * machine.beta*c/4. # RMS bunch length in meters

		self.slicer = UniformBinSlicer(
	        50, z_cuts=(-3*self.sigma_z, 3*self.sigma_z),
	        circumference=machine.circumference,
	        h_bunch=machine.h_bunch)

		# Kicker-to-pickup phase advance 0 deg
		self.damper = OneboxFeedback(self.gain,self.slicer, self.processors_detailed_x,
		                        self.processors_detailed_y, pickup_axis='displacement',
		                        kicker_axis='divergence', mpi=False,
		                        beta_x=machine.beta_x, beta_y=machine.beta_y)