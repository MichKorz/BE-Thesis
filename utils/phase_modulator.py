import sys
import json
import math
import numpy as np

C_NM_FS = 299.792458

def load_speck(filepath):
    data = np.loadtxt(filepath)
    wavelengths_nm = data[:, 0]
    
    omega = np.zeros_like(wavelengths_nm)
    valid_mask = wavelengths_nm != 0
    omega[valid_mask] = 2 * np.pi * C_NM_FS / wavelengths_nm[valid_mask]
    
    return data, wavelengths_nm, omega, data[:, 3], data[:, 4]

def apply_taylor_phase_lwe_style(omega, re, im, coefficients):
    intensity = re**2 + im**2
    
    pos_mask = omega > 0
    if not np.any(pos_mask):
        return re, im, 0.0
        
    omega_pos = omega[pos_mask]
    intensity_pos = intensity[pos_mask]
    omega_0 = omega_pos[np.argmax(intensity_pos)]
    
    d_omega = np.abs(omega) - omega_0
    abs_phase = np.zeros_like(omega)
    
    # Obliczanie wielomianu Taylora. Pętla dynamicznie skaluje się do rzędu n.
    for i, coef in enumerate(coefficients):
        order = i + 2
        abs_phase += (coef / math.factorial(order)) * (d_omega ** order)
        
    # Inwersja znaku (konwencja propagacji LWE) i zachowanie symetrii hermitowskiej
    delta_phi = -abs_phase * np.sign(omega)
    
    max_phase_shift = np.max(np.abs(delta_phi))
    
    cos_phi = np.cos(delta_phi)
    sin_phi = np.sin(delta_phi)
    
    new_re = re * cos_phi - im * sin_phi
    new_im = re * sin_phi + im * cos_phi
    
    return new_re, new_im, max_phase_shift

def write_speck(filepath, original_data, new_re, new_im):
    output_data = np.copy(original_data)
    output_data[:, 3] = new_re
    output_data[:, 4] = new_im
    np.savetxt(filepath, output_data, fmt='%1.6e')

if __name__ == "__main__":
    try:
        input_data = json.loads(sys.stdin.read())
        input_path = input_data['input_path']
        output_path = input_data['output_path']
        coefficients = input_data['coefficients']
        
        orig_data, wls_nm, omega, o_re, o_im = load_speck(input_path)
        
        new_re, new_im, max_shift = apply_taylor_phase_lwe_style(omega, o_re, o_im, coefficients)
        write_speck(output_path, orig_data, new_re, new_im)
        
        print(json.dumps({"status": "success", "max_phase_shift": max_shift}))
        sys.exit(0)
        
    except Exception as e:
        print(json.dumps({"status": "error", "message": str(e)}))
        sys.exit(1)