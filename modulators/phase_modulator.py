import numpy as np
import math

C_NM_FS = 299.792458

class PhaseModulatorBase:
    def __init__(self, input_path, output_path, **kwargs):
        self.input_path = input_path
        self.output_path = output_path

    def _load_speck(self):
        data = np.loadtxt(self.input_path)
        wavelengths_nm = data[:, 0]
        
        omega = np.zeros_like(wavelengths_nm)
        valid_mask = wavelengths_nm != 0
        omega[valid_mask] = 2 * np.pi * C_NM_FS / wavelengths_nm[valid_mask]
        
        return data, omega, data[:, 3], data[:, 4]

    def _write_speck(self, original_data, new_re, new_im):
        output_data = np.copy(original_data)
        output_data[:, 3] = new_re
        output_data[:, 4] = new_im
        np.savetxt(self.output_path, output_data, fmt='%1.6e')

    def apply_phase(self, vector):
        raise NotImplementedError()

    def get_bounds(self):
        raise NotImplementedError()


class SLMModulator(PhaseModulatorBase):
    def __init__(self, input_path, output_path, lambda_min, lambda_max, num_stripes, **kwargs):
        super().__init__(input_path, output_path, **kwargs)
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.num_stripes = num_stripes
        self.bounds = [[0.0, 2 * np.pi] for _ in range(num_stripes)]

    def apply_phase(self, vector):
        orig_data, omega, re, im = self._load_speck()
        wavelengths_nm = orig_data[:, 0]
        
        in_bounds = (wavelengths_nm >= self.lambda_min) & (wavelengths_nm <= self.lambda_max) & (wavelengths_nm != 0)
        
        new_re = np.copy(re)
        new_im = np.copy(im)
        final_phase = np.zeros_like(omega)
        
        if np.any(in_bounds):
            bins = np.linspace(self.lambda_min, self.lambda_max, self.num_stripes + 1)
            
            for i in range(self.num_stripes):
                mask = (wavelengths_nm >= bins[i]) & (wavelengths_nm <= bins[i+1]) & in_bounds
                if np.any(mask):
                    final_phase[mask] = vector[i]
                    
            cos_phi = np.cos(final_phase[in_bounds])
            sin_phi = np.sin(final_phase[in_bounds])
            
            new_re[in_bounds] = re[in_bounds] * cos_phi - im[in_bounds] * sin_phi
            new_im[in_bounds] = re[in_bounds] * sin_phi + im[in_bounds] * cos_phi
        
        self._write_speck(orig_data, new_re, new_im)
        
        orig_intensity = re**2 + im**2
        
        return wavelengths_nm, orig_intensity, final_phase

    def get_bounds(self):
        return self.bounds


class TaylorSLMModulator(PhaseModulatorBase):
    def __init__(self, input_path, output_path, lambda_min, lambda_max, num_stripes, lambda_0, bounds, **kwargs):
        super().__init__(input_path, output_path, **kwargs)
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.num_stripes = num_stripes
        self.lambda_0 = lambda_0
        self.bounds = bounds

    def apply_phase(self, vector):
        orig_data, omega, re, im = self._load_speck()
        wavelengths_nm = orig_data[:, 0]
        
        omega_0 = 2 * np.pi * C_NM_FS / self.lambda_0
        
        in_bounds = (wavelengths_nm >= self.lambda_min) & (wavelengths_nm <= self.lambda_max) & (wavelengths_nm != 0)
        
        new_re = np.copy(re)
        new_im = np.copy(im)
        final_phase = np.zeros_like(omega)
        
        if np.any(in_bounds):
            bins = np.linspace(self.lambda_min, self.lambda_max, self.num_stripes + 1)
            
            for i in range(self.num_stripes):
                mask = (wavelengths_nm >= bins[i]) & (wavelengths_nm <= bins[i+1]) & in_bounds
                if np.any(mask):
                    lambda_center = (bins[i] + bins[i+1]) / 2.0
                    omega_center = 2 * np.pi * C_NM_FS / lambda_center
                    
                    delta_omega = omega_center - omega_0
                    taylor_phase = 0.0
                    
                    for order_idx, coeff in enumerate(vector):
                        n = order_idx + 2 
                        taylor_phase += (coeff / math.factorial(n)) * (delta_omega ** n)
                        
                    wrapped_phase = taylor_phase % (2 * np.pi)
                    
                    final_phase[mask] = wrapped_phase
                    
            cos_phi = np.cos(final_phase[in_bounds])
            sin_phi = np.sin(final_phase[in_bounds])
            
            new_re[in_bounds] = re[in_bounds] * cos_phi - im[in_bounds] * sin_phi
            new_im[in_bounds] = re[in_bounds] * sin_phi + im[in_bounds] * cos_phi
        
        self._write_speck(orig_data, new_re, new_im)
        
        orig_intensity = re**2 + im**2
        
        return wavelengths_nm, orig_intensity, final_phase

    def get_bounds(self):
        return self.bounds