#!/usr/bin/env python3
"""
Copyright 2025 Thomas Berglund (NRK)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---
HLG Display Adaptation Calculator

Version: 1.3
Last modified: 2025-04-28
Author: Thomas Berglund, with assistance from Claude 3.7 Sonnet, by Anthropic.

This HLG Display Adaptation Calculator implements the Hybrid Log-Gamma (HLG) system as specified in
ITU-R BT.2100, with calculations for System Gamma and HDR Reference White following ITU-R BT.2100 Note 5f,
ITU-R BT.2390 Section 6.2 and ITU-R BT.2408 Section 2.1 respectively.

The intention is to help users understand more about how the Hybrid Log-Gamma (HLG) system works.

The idea for this tool was inspired by work related to:
    - EBU Recommendation R 167 — Reference Monitors: Predefined Modes for HLG
        https://tech.ebu.ch/publications/r167

Calculation Modes:
    nominal  : Calculates System Gamma and HDR Reference White for a standard display setup
    extended : Calculates nominal luminance of 100% and the extended 109% luminance value

Required Parameters:
    --peak     : Display peak luminance in cd/m² (default: 1000)
    --surround : Surround luminance in cd/m² (default: 5)
    --mode     : Calculation mode (default: nominal)

Optional Parameters:
    --black   : Enable black level lift with specified value in cd/m² (default: 0.0)
    --pluge   : Calculate PLUGE test signal values according to ITU-R BT.814
    --explain : Show detailed explanations of calculation steps

Usage examples:
    python3 hlg_calculator.py --peak 1000 --surround 5 --mode nominal --explain
    python3 hlg_calculator.py --peak 1000 --surround 5 --pluge --explain
    python3 hlg_calculator.py --peak 1000 --surround 5 --black 0.005 --pluge --mode nominal --explain
    python3 hlg_calculator.py --peak 1000 --surround 5 --mode extended --explain

References:
    - ITU-R BT.2100-3: Image parameter values for HDR television
        https://www.itu.int/rec/R-REC-BT.2100

    - ITU-R BT.2390: High dynamic range television for production and international programme exchange
        https://www.itu.int/pub/R-REP-BT.2390

    - ITU-R BT.2408: Operational practices in HDR television production
        https://www.itu.int/pub/R-REP-BT.2408

    - ITU-R BT.814: 
        Specifications of PLUGE test signals and alignment procedures
        for setting of brightness and contrast of displays
        https://www.itu.int/rec/R-REC-BT.814

    - High Dynamic Range Television and Hybrid Log-Gamma
        https://www.bbc.co.uk/rd/projects/high-dynamic-range

    - WHP 283 (2014): Non-linear OETF for HDR television
        https://www.bbc.co.uk/rd/publications/whitepaper283

    - WHP 309 (2015): Display-independent HDR television system
        https://www.bbc.co.uk/rd/publications/whitepaper309

    - WHP 341 (2019): HDR TV brightness perception modeling
        https://www.bbc.co.uk/rd/publications/whitepaper341

    - WHP 342 (2019): HDR brightness perception and transitions
        https://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP342.pdf

    - WHP 369 (2019):Display of High Dynamic Range Images Under Varying Viewing Conditions
        https://www.bbc.co.uk/rd/publications/display-high-dynamic-range-images-varying-viewing-conditions
"""

import math
import argparse

# Global flag to control explanation output
# (default: minimal output with only calculated results)
EXPLAIN = False

def e_print(*args, **kwargs):
    """Print detailed explanation steps only when the explain flag is active."""
    if EXPLAIN:
        print(*args, **kwargs)

"""
ITU-R BT.2100-3 HLG Constants

These constants are defined explicitly according to ITU-R BT.2100-3 Table 5.
They are used in the transfer function calculations for the HLG system.
"""
# HLG constants defined in ITU-R BT.2100-3 Table 5
HLG_a = 0.17883277
HLG_b = 1 - 4 * HLG_a  # Equals 0.28466892
HLG_c = 0.5 - HLG_a * math.log(4 * HLG_a)  # Equals 0.55991073

# BT.2020 luminance coefficients as defined in ITU-R BT.2100-3 Table 5
Y_R = 0.2627  # Red luminance coefficient
Y_G = 0.6780  # Green luminance coefficient 
Y_B = 0.0593  # Blue luminance coefficient

# ITU-R BT.2100-3 defined 109% signal level (1.090182648401826)
EXTENDED_RANGE_109 = (1019 - 64) / (940 - 64)

#-----------------------------------------------------------------------------
# Section 1: HLG OETF and Inverse OETF
#-----------------------------------------------------------------------------

def hlg_reference_oetf(E):
    """
    HLG Opto-Electronic Transfer Function (OETF) as defined in ITU-R BT.2100-3 Table 5.
    
    Transforms scene-linear light to the non-linear HLG signal.
    
    Args:
        E (float): Scene linear light normalized to the range [0:1]
        
    Returns:
        float: Non-linear HLG signal value E' in the range [0:1]
        
    Reference: ITU-R BT.2100-3, Table 5, HLG Reference OETF
    """
    # Implement the piecewise HLG OETF function
    if E <= 1/12:
        return math.sqrt(3 * E)
    else:
        return HLG_a * math.log(12 * E - HLG_b) + HLG_c

def hlg_inverse_oetf(E_prime, explain_label=None):
    """
    HLG Inverse OETF: Converts non-linear signal value to scene linear light
    as defined in ITU-R BT.2100-3 Table 5.
    
    This is the first step of the EOTF, converting the encoded signal back to
    scene-referred linear light before applying the OOTF.
    
    Args:
        E_prime (float): Non-linear HLG signal value (E'), clipped to ≥ 0.
                         Values typically in range [0:1] but can exceed 1.0
                         in the super-white region per BT.2100-3 Note 5h.
        explain_label (str, optional): Label for explanation output
    
    Returns:
        float: Scene linear light (E) normalized to the range [0:1].
               A value of 1/12 corresponds to the transition point (E'=0.5).
    
    Reference: ITU-R BT.2100-3, Table 5, Hybrid Log-Gamma (HLG) system,
               "inverse OETF" component of the "Reference HLG EOTF" definition.
    """
    # Ensure signal is non-negative per BT.2100-3 Table 5
    E_prime = max(0.0, E_prime)
    
    # Apply the piecewise inverse OETF
    if E_prime < 0.5:
        # Square law segment (E ≤ 1/12 in original OETF)
        E = (E_prime ** 2) / 3.0
        
        # Explain the calculation if requested
        if explain_label:
            e_print(f"    Applying inverse OETF to convert signal to scene linear light:")
            e_print(f"      Non-linear signal value: {E_prime:.6f}")
            e_print(f"      Signal < 0.5: Using square law formula for inverse OETF")
            e_print(f"      Formula: E = (E')² / 3")
            e_print(f"      E = ({E_prime:.6f})² / 3 = {E_prime**2:.6f} / 3 = {E:.6f}")
    else:
        # Logarithmic segment (E > 1/12 in original OETF)
        # Calculate intermediate values for explanation
        exp_term = math.exp((E_prime - HLG_c) / HLG_a)
        E = (exp_term + HLG_b) / 12.0
        
        # Explain the calculation if requested
        if explain_label:
            e_print(f"    Applying inverse OETF to convert signal to scene linear light:")
            e_print(f"      Non-linear signal value: {E_prime:.6f}")
            e_print(f"      Signal ≥ 0.5: Using logarithmic formula for inverse OETF")
            e_print(f"      Formula: E = (exp((E' - c) / a) + b) / 12")
            e_print(f"      Where: a = {HLG_a:.8f}, b = {HLG_b:.8f}, c = {HLG_c:.8f}")
            e_print(f"      E = (exp(({E_prime:.6f} - {HLG_c:.8f}) / {HLG_a:.8f}) + {HLG_b:.8f}) / 12")
            e_print(f"      E = (exp({(E_prime - HLG_c):.8f} / {HLG_a:.8f}) + {HLG_b:.8f}) / 12")
            e_print(f"      E = (exp({(E_prime - HLG_c) / HLG_a:.8f}) + {HLG_b:.8f}) / 12")
            e_print(f"      E = ({exp_term:.8f} + {HLG_b:.8f}) / 12 = {E:.6f}")
    
    return E

#-----------------------------------------------------------------------------
# Section 2: HLG System Gamma and Black Level Lift
#-----------------------------------------------------------------------------

def calculate_system_gamma(Lw, Ls, step_number=None):
    """
    Calculate the System Gamma (γ) according to formulas specified in:
    - ITU-R BT.2100-3 Note 5f
    - ITU-R BT.2390 Section 6.2
    
    The System Gamma is adapted according to the two following conditions:
    1. The display's peak luminance capability (Lw)
    2. The viewing environment's surround luminance (Ls)
    
    Two different formulas are specified in ITU-R BT.2100-3 Note 5f.
    The recommended formula depends on the display's peak luminance,
    each optimized for different peak display luminance ranges:

    - For 400 cd/m² ≤ Lw ≤ 2000 cd/m²: γ = 1.2 + 0.42 log₁₀(Lw/1000)
    - For Lw < 400 cd/m² or Lw > 2000 cd/m²: γ = 1.2 × κ^(log₂(Lw/1000))
    
    ITU-R BT.2100-3 Note 5f:
    "For displays with nominal peak luminance (Lw) other than 1 000 cd/m2,
    or where the effective nominal peak luminance is adjusted through the use
    of a contrast control, the System Gamma value should be adjusted according
    to the formula below and may be rounded to three significant digits:"
    
        γ = 1.2 + 0.42 log₁₀(Lw/1000)

    "For applications outside of the usual production monitoring range of Lw equal to
    400 cd/m² to 2 000 cd/m², the following extended range formula should be used:"

        γ = 1.2 × κ^(log₂(Lw/1000)) where κ = 1.111
        
    The formula is centered around the reference point of 1000 cd/m²,
    where both formulas produce identical results with a System Gamma of γ = 1.2
    
    In both cases, additional adjustment for surround luminance should be applied:
        γ = γ × μ^(log₂(Ls/1000))
        
    The complete extended formula from ITU-R BT.2390 is:
        γ = γ_ref × κ^(log₂(Lw/L_ref)) × μ^(log₂(Ls/Ls_ref))
        
    Where:
        - γ_ref = 1.2 (reference gamma at 1000 cd/m² in reference environment)
        - κ = 1.111 (peak luminance adjustment factor)
        - μ = 0.98 (surround luminance adjustment factor)
        - L_ref = 1000 cd/m² (reference display peak luminance)
        - Ls_ref = 5 cd/m² (reference surround luminance)

    Args:
        Lw (float): Nominal peak luminance of the display in cd/m².
        Ls (float): Surround luminance in cd/m².
        step_number (int, optional): Step number for explanation purposes.
    
    Returns:
        tuple: (system_gamma, calculation_details_dict)
              Where calculation_details_dict contains all intermediate values
    
    Note:
        Higher display peak luminance requires increased gamma (κ > 1)
        Brighter viewing environments require decreased gamma (μ < 1)
    
    Reference: ITU-R BT.2100-3 Note 5f and ITU-R BT.2390 Section 6.2
    """
    # Constants
    gamma_ref = 1.2  # γ_ref
    kappa = 1.111    # κ
    mu = 0.98        # μ
    L_ref = 1000     # L_ref
    Ls_ref = 5       # Ls_ref
    
    # Select the appropriate formula based on display peak luminance
    if 400 <= Lw <= 2000:
        # Basic formula for typical production range (400-2000 cd/m²)
        # γ = 1.2 + 0.42 log₁₀(Lw/L_ref)
        base_gamma = gamma_ref + 0.42 * math.log10(Lw / L_ref)
        kappa_factor = None  # Not used with basic formula
    else:
        # Extended formula for displays outside typical range
        # γ = 1.2 × κ^(log₂(Lw/L_ref))
        kappa_factor = kappa ** math.log2(Lw / L_ref)
        base_gamma = gamma_ref * kappa_factor
    
    # Calculate surround luminance adjustment factor (reduces gamma for brighter environments)
    # Following ITU-R BT.2390 Section 6.2
    mu_factor = mu ** math.log2(Ls / Ls_ref)
    
    # Apply surround adjustment to System Gamma
    gamma = base_gamma * mu_factor
    
    # Print explanation if step number is provided
    if step_number is not None:
        e_print(f"\nStep {step_number}: Calculating System Gamma (γ)")
        e_print(f"  Peak luminance (Lw): {Lw:.2f} cd/m²")
        e_print(f"  Surround luminance (Ls): {Ls:.2f} cd/m²")
        
        if 400 <= Lw <= 2000:
            e_print(f"  Using basic formula for 400 ≤ Lw ≤ 2000 cd/m²:")
            e_print(f"  Formula: γ = 1.2 + 0.42 log₁₀(Lw/1000)")
            e_print(f"  Base gamma: 1.2 + 0.42 × log₁₀({Lw:.2f}/1000) = {base_gamma:.6f}")
        else:
            e_print(f"  Using extended formula for Lw < 400 or Lw > 2000 cd/m²:")
            e_print(f"  Formula: γ = γref × κ^(log₂(Lw/L_ref))")
            e_print(f"  Where: γref = {gamma_ref}, κ = {kappa}")
            e_print(f"  Peak luminance adjustment: κ^(log₂({Lw:.2f}/{L_ref})) = {kappa_factor:.6f}")
            e_print(f"  Base gamma: {gamma_ref} × {kappa_factor:.6f} = {base_gamma:.6f}")
        e_print()
        e_print(f"  Surround adjustment formula: μ^(log₂(Ls/Ls_ref))")
        e_print(f"  Where: μ = {mu}, Ls_ref = {Ls_ref} cd/m²")
        e_print(f"  Surround luminance adjustment: μ^(log₂({Ls:.2f}/{Ls_ref})) = {mu_factor:.6f}")
        e_print(f"  Final System Gamma (γ): {base_gamma:.6f} × {mu_factor:.6f} = {gamma:.6f}")
        
        # Add note if reference surround
        if Ls == 5.0:
            e_print(f"  Note: {Ls:.2f} cd/m² is the reference surround luminance")
    
    # Return both the System Gamma and calculation details
    return gamma, {
        'gamma_ref': gamma_ref,   # γ_ref
        'kappa': kappa,           # κ
        'mu': mu,                 # μ
        'L_ref': L_ref,           # L_ref
        'Ls_ref': Ls_ref,         # Ls_ref
        'base_gamma': base_gamma,
        'kappa_factor': kappa_factor,
        'mu_factor': mu_factor
    }

def calculate_black_level_lift(L_B, L_W, gamma, step_number=None):
    """
    Calculate the black level lift (β) according to ITU-R BT.2100-3.
    
    The black level lift parameter adjusts signal mapping to accommodate non-zero
    black levels in real displays, ensuring consistent rendering across displays
    with different black level capabilities.
    
    Per ITU-R BT.2100-3 Table 5, black level lift is defined as:
        β = √3(LB/LW)^(1/γ)
    
    Where:
    - LB is the display luminance for black in cd/m²
    - LW is nominal peak luminance of the display in cd/m²
    - γ is the System Gamma
    
    This parameter ensures proper rendering of the signal's black level on displays
    where absolute zero luminance isn't achievable. It modifies the signal mapping
    to maintain perceptual consistency across different display technologies.
    
    Args:
        L_B (float): Display luminance for black in cd/m².
        L_W (float): Nominal peak luminance of the display in cd/m².
        gamma (float): System gamma calculated for the display and viewing environment.
        step_number (int, optional): Step number for explanation purposes.
        
    Returns:
        float: Black level lift factor (β).
    
    Reference: ITU-R BT.2100-3, Table 5, HLG Reference EOTF, definition of β.
    """
    # Early return if black level is 0 (ideal display)
    if L_B <= 0:
        if step_number is not None:
            e_print(f"\nStep {step_number}: Calculating Black Level Lift (β)")
            e_print(f"  Black level (L_B) is {L_B:.4f} cd/m²")
            e_print(f"  Since black level is zero or negative, black level lift is disabled (β = 0)")
        return 0.0
        
    # Calculate black level lift: β = √3(LB/LW)^(1/γ)
    beta = math.sqrt(3 * ((L_B / L_W) ** (1 / gamma)))
    
    # Explanation output
    if step_number is not None:
        e_print(f"\nStep {step_number}: Calculating Black Level Lift (β)")
        e_print(f"  Black level (L_B): {L_B:.4f} cd/m²")
        e_print(f"  Peak luminance (L_W): {L_W:.2f} cd/m²")
        e_print(f"  System gamma (γ): {gamma:.4f}")
        e_print(f"  Formula: β = √3(L_B/L_W)^(1/γ)")
        e_print(f"  β = √3({L_B:.4f}/{L_W:.2f})^(1/{gamma:.4f}) = {beta:.6f}")
        
        # Educational notes about black level lift
        e_print(f"  Note: Black level lift compensates for the display's inability")
        e_print(f"        to produce perfect black (zero luminance).")
        e_print(f"  Note: Per ITU-R BT.2100-3 Note 5h, values below 0.0 should not be clipped")
        e_print(f"        to allow proper black level setting with PLUGE test signals.")
    
    return beta

#-----------------------------------------------------------------------------
# Section 3: HLG Reference OOTF
#-----------------------------------------------------------------------------

def hlg_reference_ootf(r_s, g_s, b_s, gamma, alpha=1.0, explain_label=None):
    """
    HLG Reference Opto-Optical Transfer Function (OOTF) as defined in ITU-R BT.2100-3 Table 5.
    
    Maps scene-linear RGB to display-linear RGB with appropriate gamma adjustment.
    The HLG OOTF applies the System Gamma to the luminance component only,
    preserving consistent color appearance through the rendering chain.
    
    As specified in ITU-R BT.2390 Section 6.2:
    "Instead of the current SDR practice of applying a gamma curve independently
    to each colour component, for HDR it should be applied to the luminance alone."
    
    The HLG OOTF calculates a single luminance adjustment factor (Y_S^(γ-1))
    and applies it uniformly to all three color components. This preserves the
    ratios between R, G, and B, maintaining consistent color appearance across
    different luminance levels.
    
    Args:
        r_s, g_s, b_s (float): Scene-linear RGB components [0:1]
        gamma (float): System gamma calculated based on display and environment
        alpha (float, optional): Display nominal peak luminance scaling factor.
                                 Typically set to 1.0 for normalized calculations.
        explain_label (str, optional): Label for explanation output
    
    Returns:
        tuple: Display-linear RGB components (r_d, g_d, b_d)
    
    Reference: ITU-R BT.2100-3, Table 5, HLG Reference OOTF
    """
    # Calculate scene luminance Y_S using BT.2020 coefficients
    # ITU-R BT.2100-3 Table 5: Y_S = 0.2627R_S + 0.6780G_S + 0.0593B_S
    Y_S = Y_R * r_s + Y_G * g_s + Y_B * b_s
    
    # Calculate OOTF luminance adjustment factor (Y_S^(γ-1))
    # This is the critical step where gamma is applied to luminance, not to individual components
    if Y_S > 0:
        ootf_factor = Y_S ** (gamma - 1)
    else:
        ootf_factor = 0  # Avoid math domain error for negative/zero values
    
    # Apply OOTF luminance adjustment to each component
    # This preserves color ratios (chromaticity) while adjusting overall luminance
    r_d = alpha * ootf_factor * r_s  # Display-referred red
    g_d = alpha * ootf_factor * g_s  # Display-referred green
    b_d = alpha * ootf_factor * b_s  # Display-referred blue
    
    # Generate detailed explanation of HLG OOTF calculations for achromatic signals
    # where the math simplifies and becomes more intuitive to understand
    if explain_label and r_s == g_s == b_s:
        # For achromatic signals (R=G=B), calculations are more straightforward
        # and easier to follow, making them ideal for educational purposes
        L = r_s  # Scene linear light value (identical for all components)
        
        e_print(f"  For achromatic signal (R=G=B, {explain_label}):")
        
        # Explain scene luminance calculation
        e_print(f"    Scene luminance Y_S = {Y_R:.4f}*{r_s:.4f} + {Y_G:.4f}*{g_s:.4f} + {Y_B:.4f}*{b_s:.4f} = {Y_S:.4f}")
        
        # Show calculation of the luminance adjustment factor per ITU-R BT.2100-3 Table 5
        e_print(f"    OOTF luminance adjustment factor = Y_S^(γ-1) = {Y_S:.4f}^{gamma - 1:.4f} = {ootf_factor:.4f}")
        
        # Demonstrate how the OOTF preserves the same adjustment across all three components
        e_print(f"    Display light F_D = OOTF[L] = α × L × Y_S^(γ-1) = {alpha:.4f} × {L:.4f} × {ootf_factor:.4f} = {r_d:.4f}")
        
        # Show the mathematical simplification that occurs with achromatic signals
        e_print(f"    Note: For achromatic signals, this simplifies to F_D = α × L^γ = {alpha:.4f} × {L:.4f}^{gamma:.4f} = {r_d:.4f}")
    
    return r_d, g_d, b_d

#-----------------------------------------------------------------------------
# Section 4: HLG Reference EOTF (using Inverse OETF + OOTF)
#-----------------------------------------------------------------------------

def hlg_reference_eotf(r_prime, g_prime, b_prime, gamma, L_W=None, L_B=0.0, beta=None, explain_label=None):
    """
    Implements the full HLG Electro-Optical Transfer Function (EOTF) with OOTF 
    for RGB signals as specified in ITU-R BT.2100-3 Table 5, with optional
    black level lift support.
    
    The HLG Reference EOTF consists of these conceptual steps:
    1. Apply black level lift if required: E' → max(0,(1-β)E'+β)
    2. Inverse OETF: Convert non-linear R', G', B' to scene linear R, G, B
    3. Apply OOTF: F_D = α·Y_S^(γ-1)·E, where:
       - Y_S is the scene luminance: 0.2627R + 0.6780G + 0.0593B
       - γ is the variable System Gamma
       - α is the display nominal peak luminance
    
    The complete chain is expressed as: F_D = OOTF[OETF⁻¹[max(0,(1-β)E'+β)]]
    
    Args:
        r_prime, g_prime, b_prime (float): Non-linear HLG signal values [0:1+]
        gamma (float): System gamma calculated according to ITU-R BT.2100-3 Note 5f
                       Typically in range 1.0-1.5, based on display peak luminance
        L_W (float, optional): Nominal peak luminance of the display in cd/m²
                              Required for black level lift calculation
        L_B (float, optional): Display luminance for black in cd/m²
        beta (float, optional): Precalculated black level lift factor to reuse
                               If None, will be calculated from L_B and L_W
        explain_label (str, optional): Label for explanation output
    
    Returns:
        tuple: Normalized display light values (r_d, g_d, b_d)
               Multiply by nominal peak luminance to get actual display luminance
    
    Reference: ITU-R BT.2100-3, Table 5, "HLG Reference EOTF" definition.
    """
    # Step 1: Apply black level lift if required (L_B > 0)
    black_level_enabled = (L_B > 0.0 and L_W is not None) or (beta is not None and beta > 0)
    
    if black_level_enabled:
        # Use provided beta or calculate it
        if beta is None:
            beta = calculate_black_level_lift(L_B, L_W, gamma)
        
        # Apply black level lift to each component per BT.2100-3 formula:
        # F_D = EOTF[max(0,(1-β)E'+β)]
        r_lifted = max(0, (1 - beta) * r_prime + beta)
        g_lifted = max(0, (1 - beta) * g_prime + beta)
        b_lifted = max(0, (1 - beta) * b_prime + beta)
        
        # Store original values for explanation
        if explain_label and r_prime == g_prime == b_prime:
            e_print(f"  Applying black level lift (β = {beta:.6f}):")
            e_print(f"    Original signal: {r_prime:.6f}")
            e_print(f"    Formula: E'_lifted = max(0, (1-β)E' + β)")
            e_print(f"    Lifted signal: max(0, (1-{beta:.6f})*{r_prime:.6f} + {beta:.6f}) = {r_lifted:.6f}")
    else:
        # No black level lift
        r_lifted = r_prime
        g_lifted = g_prime
        b_lifted = b_prime
    
    # Step 2: Convert non-linear signal to scene linear light using inverse OETF
    # Only pass the explain_label for achromatic signals to avoid redundant explanations
    explain_inverse_oetf = explain_label and r_lifted == g_lifted == b_lifted
    
    r_s = hlg_inverse_oetf(r_lifted, explain_label if explain_inverse_oetf else None)  # Scene linear red
    g_s = hlg_inverse_oetf(g_lifted)  # Scene linear green
    b_s = hlg_inverse_oetf(b_lifted)  # Scene linear blue
    
    # Step 3: Apply the OOTF to get display-referred light
    # Using the separated OOTF function (α=1.0 for normalized values)
    r_d, g_d, b_d = hlg_reference_ootf(r_s, g_s, b_s, gamma, alpha=1.0, explain_label=explain_label)
    
    return r_d, g_d, b_d

#-----------------------------------------------------------------------------
# Section 5: HDR Reference White and PLUGE
#-----------------------------------------------------------------------------

def calculate_reference_white(Lw, gamma, L_B=0.0, beta=None, step_number=None):
    """
    Calculate the luminance of HDR Reference White (75% signal level).

    Based on ITU-R BT.2408 Section 2.1, HDR Reference White is defined as:
    "the nominal signal level obtained from an HDR camera and a 100% reflectance 
    white card resulting in a nominal luminance of 203 cd/m² on a PQ display or 
    on an HLG display that has a nominal peak luminance capability of 1,000 cd/m²."
    
    For a reference 1000 cd/m² display with standard gamma of 1.2, this function 
    will return 203 cd/m², exactly matching the reference level in ITU-R BT.2408 Table 1.
    
    With black level lift enabled, the function accounts for the modified
    signal mapping as specified in ITU-R BT.2100-3.
    
    Args:
        Lw (float): Nominal peak luminance of the display in cd/m².
        gamma (float): System gamma calculated for the display and viewing environment.
        L_B (float, optional): Display luminance for black in cd/m².
        beta (float, optional): Precalculated black level lift to reuse.
        step_number (int, optional): Step number for explanation purposes.
    
    Returns:
        float: Luminance of HDR Reference White in cd/m².
    
    Reference: ITU-R BT.2408 Section 2.1
    """
    signal_75 = 0.75
    
    if step_number is not None:
        e_print(f"\nStep {step_number}: Calculating HDR Reference White (75% HLG Signal)")
    
    e_print(f"  Input signal value: {signal_75:.2f} (75%)")
    
    # Calculate normalized luminance using the HLG EOTF with black level lift if enabled
    r_d, _, _ = hlg_reference_eotf(signal_75, signal_75, signal_75, gamma, 
                         L_W=Lw, L_B=L_B, beta=beta, explain_label="75%")
    L_normalized = r_d
    
    e_print(f"  Normalized luminance: {L_normalized:.6f}")
    
    ref_white = Lw * L_normalized
    e_print(f"  Multiply by nominal peak: {Lw:.6f} × {L_normalized:.6f} = {ref_white:.6f} cd/m²")
    
    return ref_white

def calculate_pluge_values(peak_luminance, surround_luminance, black_level=0.0):
    """
    Calculate PLUGE test signal levels for HLG displays according to ITU-R BT.814.
    
    This function calculates the luminance values corresponding to the 10-bit and 12-bit
    code values specified in ITU-R BT.814 Table 3 for HDR displays. It uses the HLG EOTF
    with System Gamma adapted to the viewing environment.
    
    Args:
        peak_luminance (float): Nominal peak luminance of the display in cd/m².
        surround_luminance (float): Surround luminance in cd/m².
        black_level (float, optional): Display luminance for black in cd/m².
        
    Returns:
        dict: A dictionary containing PLUGE level luminances and normalization factors:
            - 'higher_level': Luminance of the Higher level (code value 399/1596)
            - 'black_level_display': Luminance of the Black level (code value 64/256)
            - 'slightly_lighter_level': Luminance of Slightly lighter level (code value 80/320)
            - 'slightly_darker_level': Luminance of Slightly darker level (code value 48/192)
            - 'system_gamma': The calculated System Gamma
            - 'black_level_lift': The calculated black level lift factor (β)
            - 'input_black_level': Original black level input parameter
    """
    # Store all values to return instead of printing during calculation
    gamma = None
    beta = None
    higher_level_luminance = None
    black_level_display_luminance = None
    slightly_lighter_luminance = None
    slightly_darker_luminance = None
    
    # Only print header if in explain mode
    if EXPLAIN:
        e_print("\nPLUGE Test Signal Calculation:")
        e_print(f"{'-' * 75}")
        e_print(f"Input Peak Luminance (Lw): {peak_luminance:.2f} cd/m²")
        e_print(f"Input Surround Luminance (Ls): {surround_luminance:.2f} cd/m²")
        e_print(f"Input Black Level (LB): {black_level:.4f} cd/m²")
    
    # Step 1: Calculate System Gamma
    gamma, _ = calculate_system_gamma(peak_luminance, surround_luminance, step_number=1 if EXPLAIN else None)
    
    # Step 2: Calculate Black Level Lift (if enabled)
    beta = 0.0
    if black_level > 0.0:
        beta = calculate_black_level_lift(black_level, peak_luminance, gamma, step_number=2 if EXPLAIN else None)
    
    # Define code values from ITU-R BT.814 Table 3
    # 10-bit values
    higher_level_code_10bit = 399
    black_level_code_10bit = 64
    slightly_lighter_code_10bit = 80
    slightly_darker_code_10bit = 48
    
    # 12-bit values
    higher_level_code_12bit = 1596
    black_level_code_12bit = 256
    slightly_lighter_code_12bit = 320
    slightly_darker_code_12bit = 192
    
    # Define normalization parameters for different bit depths
    # 10-bit normalization parameters
    black_10bit = 64
    peak_10bit = 940
    
    # 12-bit normalization parameters
    black_12bit = 256
    peak_12bit = 3760
    
    # Step 3: Convert code values to normalized signal values
    if EXPLAIN:
        e_print(f"\nStep {3 if black_level <= 0 else 3}: Converting code values to normalized signal values")
    
    # Calculate normalized values from 10-bit codes
    higher_level_normalized = (higher_level_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    black_level_normalized = (black_level_code_10bit - black_10bit) / (peak_10bit - black_10bit)  # Should be 0
    slightly_lighter_normalized = (slightly_lighter_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    
    # For slightly darker (which is below black level), calculate the negative offset
    slightly_darker_normalized = (slightly_darker_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    
    if EXPLAIN:
        # Show normalization for 10-bit codes
        e_print(f"  10-bit Higher level code: {higher_level_code_10bit} → normalized: {higher_level_normalized:.4f}")
        e_print(f"  10-bit Black level code: {black_level_code_10bit} → normalized: {black_level_normalized:.4f}")
        e_print(f"  10-bit Slightly lighter code: {slightly_lighter_code_10bit} → normalized: {slightly_lighter_normalized:.4f}")
        e_print(f"  10-bit Slightly darker code: {slightly_darker_code_10bit} → normalized: {slightly_darker_normalized:.4f}")
        
        # Verify with 12-bit codes (should give same normalized values)
        higher_level_normalized_12bit = (higher_level_code_12bit - black_12bit) / (peak_12bit - black_12bit)
        black_level_normalized_12bit = (black_level_code_12bit - black_12bit) / (peak_12bit - black_12bit)
        slightly_lighter_normalized_12bit = (slightly_lighter_code_12bit - black_12bit) / (peak_12bit - black_12bit)
        slightly_darker_normalized_12bit = (slightly_darker_code_12bit - black_12bit) / (peak_12bit - black_12bit)
        
        # Show normalization for 12-bit codes
        e_print(f"  12-bit Higher level code: {higher_level_code_12bit} → normalized: {higher_level_normalized_12bit:.4f}")
        e_print(f"  12-bit Black level code: {black_level_code_12bit} → normalized: {black_level_normalized_12bit:.4f}")
        e_print(f"  12-bit Slightly lighter code: {slightly_lighter_code_12bit} → normalized: {slightly_lighter_normalized_12bit:.4f}")
        e_print(f"  12-bit Slightly darker code: {slightly_darker_code_12bit} → normalized: {slightly_darker_normalized_12bit:.4f}")
        
        # Verify both normalizations match
        e_print(f"  Both 10-bit and 12-bit normalization should match: {higher_level_normalized:.4f} ≈ {higher_level_normalized_12bit:.4f}")
    
    # Step 4: Calculate luminance values using the HLG EOTF
    if EXPLAIN:
        e_print(f"\nStep {4 if black_level <= 0 else 4}: Calculating PLUGE Level Luminances")
    
    # Calculate Higher level luminance (the reference 38.2% level)
    if EXPLAIN:
        e_print(f"  Higher level (code value 399/1596, normalized signal {higher_level_normalized:.3f}):")
    r_d, _, _ = hlg_reference_eotf(
        higher_level_normalized, 
        higher_level_normalized, 
        higher_level_normalized, 
        gamma, 
        L_W=peak_luminance, 
        L_B=black_level, 
        beta=beta,
        explain_label="Higher level" if EXPLAIN else None
    )
    higher_level_luminance = peak_luminance * r_d
    if EXPLAIN:
        e_print(f"  Higher level luminance = peak_luminance × normalized_value = {peak_luminance:.2f} × {r_d:.4f} = {higher_level_luminance:.2f} cd/m²")
    
    # Calculate Black level luminance
    if EXPLAIN:
        e_print(f"  Black level (code value 64/256, normalized signal {black_level_normalized:.3f}):")
    r_d, _, _ = hlg_reference_eotf(
        black_level_normalized, 
        black_level_normalized, 
        black_level_normalized, 
        gamma, 
        L_W=peak_luminance, 
        L_B=black_level, 
        beta=beta,
        explain_label="Black level" if EXPLAIN else None
    )
    black_level_display_luminance = peak_luminance * r_d
    if EXPLAIN:
        e_print(f"  Black level luminance = peak_luminance × normalized_value = {peak_luminance:.2f} × {r_d:.4f} = {black_level_display_luminance:.4f} cd/m²")
    
    # Calculate Slightly lighter level luminance
    if EXPLAIN:
        e_print(f"  Slightly lighter level (code value 80/320, normalized signal {slightly_lighter_normalized:.3f}):")
    r_d, _, _ = hlg_reference_eotf(
        slightly_lighter_normalized, 
        slightly_lighter_normalized, 
        slightly_lighter_normalized, 
        gamma, 
        L_W=peak_luminance, 
        L_B=black_level, 
        beta=beta,
        explain_label="Slightly lighter level" if EXPLAIN else None
    )
    slightly_lighter_luminance = peak_luminance * r_d
    if EXPLAIN:
        e_print(f"  Slightly lighter level luminance = peak_luminance × normalized_value = {peak_luminance:.2f} × {r_d:.4f} = {slightly_lighter_luminance:.4f} cd/m²")
    
    # Calculate Slightly darker level luminance
    if EXPLAIN:
        e_print(f"  Slightly darker level (code value 48/192, normalized signal {slightly_darker_normalized:.2f}):")
    r_d, _, _ = hlg_reference_eotf(
        slightly_darker_normalized, 
        slightly_darker_normalized, 
        slightly_darker_normalized, 
        gamma, 
        L_W=peak_luminance, 
        L_B=black_level, 
        beta=beta,
        explain_label="Slightly darker level" if EXPLAIN else None
    )
    slightly_darker_luminance = peak_luminance * r_d
    if EXPLAIN:
        e_print(f"  Slightly darker level luminance = peak_luminance × normalized_value = {peak_luminance:.2f} × {r_d:.4f} = {slightly_darker_luminance:.4f} cd/m²")
    
    if EXPLAIN:
        e_print(f"{'-' * 75}")
    
    # Return all calculated values in a dictionary
    return {
        'higher_level': higher_level_luminance,
        'black_level_display': black_level_display_luminance,
        'slightly_lighter_level': slightly_lighter_luminance,
        'slightly_darker_level': slightly_darker_luminance,
        'system_gamma': gamma,
        'black_level_lift': beta,
        'input_black_level': black_level
    }

#-----------------------------------------------------------------------------
# Section 6: Mode-specific Calculation Functions
#-----------------------------------------------------------------------------

def nominal_range(peak_luminance, surround_luminance, black_level=0.0):
    """
    Calculate standard parameters for an HLG display.
    
    This function calculates the essential parameters for an HLG display setup,
    including System Gamma adaptation based on viewing environment, optional
    black level lift, and HDR Reference White level.

    Args:
        peak_luminance (float): Nominal peak luminance of the display in cd/m².
        surround_luminance (float): Luminance of the display surround in cd/m².
        black_level (float, optional): Display luminance for black in cd/m².

    Returns:
        dict: A dictionary containing:
            - '100%_luminance' (float): Nominal 100% luminance.
            - 'system_gamma' (float): Calculated System Gamma.
            - 'black_level' (float): Display luminance for black.
            - 'black_level_lift' (float): Calculated black level lift factor (β).
            - 'reference_white' (float): HDR Reference White in cd/m².
            - 'surround_luminance' (float): Input surround luminance (for reference).
    """
    e_print("\nNominal Range Calculation Steps:")
    e_print(f"{'-' * 75}")
    e_print(f"Input Peak Luminance (Lw): {peak_luminance:.2f} cd/m²")
    e_print(f"Input Surround Luminance (Ls): {surround_luminance:.2f} cd/m²")
    e_print(f"Input Black Level (LB): {black_level:.4f} cd/m²")

    # Step 1: Calculate System Gamma using the shared function with step number
    gamma, _ = calculate_system_gamma(peak_luminance, surround_luminance, step_number=1)

    # Step 2: Calculate Black Level Lift (if enabled)
    beta = 0.0
    if black_level > 0.0:
        beta = calculate_black_level_lift(black_level, peak_luminance, gamma, step_number=2)
        # Bump the next step number
        ref_white_step = 3
    else:
        ref_white_step = 2
    
    # Step 3: Calculate HDR Reference White using the helper function with step number
    ref_white = calculate_reference_white(peak_luminance, gamma, 
                                        L_B=black_level, beta=beta, 
                                        step_number=ref_white_step)
    
    e_print(f"{'-' * 75}")
    
    return {
        '100%_luminance': peak_luminance,
        'system_gamma': gamma,
        'black_level': black_level,
        'black_level_lift': beta,
        'reference_white': ref_white,
        'surround_luminance': surround_luminance
    }

def extended_range(peak_luminance, surround_luminance, black_level=0.0):
    """
    Calculate extended range parameters for an HLG display, showing the brightness
    at 109% signal level (super-white) for a given nominal 100% luminance.
    
    This function provides detailed calculation steps to show how HLG handles signals
    beyond the nominal range, with optional support for black level lift.
    
    Args:
        peak_luminance (float): Nominal peak luminance of the display in cd/m²
        surround_luminance (float): Luminance of the display surround in cd/m²
        black_level (float, optional): Display luminance for black in cd/m²
        
    Returns:
        dict: A dictionary containing:
            - '100%_luminance' (float): Nominal 100% luminance.
            - '109%_luminance' (float): Calculated extended luminance at 109% signal.
            - 'system_gamma' (float): Calculated System Gamma.
            - 'black_level' (float): Display luminance for black.
            - 'black_level_lift' (float): Calculated black level lift factor (β).
            - 'reference_white' (float): HDR Reference White in cd/m².
    """
    e_print("\nExtended Range Calculation Steps:")
    e_print(f"{'-' * 75}")
    e_print(f"Input Peak Luminance (Lw): {peak_luminance:.2f} cd/m²")
    e_print(f"Input Surround Luminance (Ls): {surround_luminance:.2f} cd/m²")
    e_print(f"Input Black Level (LB): {black_level:.4f} cd/m²")
    
    # Step 1: Calculate System Gamma using the shared function with step number
    gamma, _ = calculate_system_gamma(peak_luminance, surround_luminance, step_number=1)
    
    # Step 2: Calculate Black Level Lift (if enabled)
    beta = 0.0
    step_offset = 0
    if black_level > 0.0:
        beta = calculate_black_level_lift(black_level, peak_luminance, gamma, step_number=2)
        step_offset = 1
    
    # Step 3: Calculate 109% signal luminance
    e_print(f"\nStep {2 + step_offset}: Calculating 109% signal luminance")
    E = EXTENDED_RANGE_109
    e_print(f"  Input signal value: {E:.15f} (109%)")
    
    # Calculate normalized luminance using HLG EOTF with black level lift if enabled
    r_d, _, _ = hlg_reference_eotf(E, E, E, gamma, L_W=peak_luminance, L_B=black_level, 
                       beta=beta, explain_label="109%")
    L_normalized = r_d
    
    e_print(f"  Normalized luminance: {L_normalized:.6f}")
    
    # Calculate luminance value
    luminance_109 = peak_luminance * L_normalized
    e_print(f"  Multiply by nominal peak: {peak_luminance:.6f} × {L_normalized:.6f} = {luminance_109:.6f} cd/m²")
    
    # Step 4: Calculate HDR Reference White using the helper function with step number
    ref_white = calculate_reference_white(peak_luminance, gamma, L_B=black_level, 
                                        beta=beta, step_number=3 + step_offset)
    
    e_print(f"{'-' * 75}")
    
    return {
        '100%_luminance': peak_luminance,
        '109%_luminance': luminance_109,
        'system_gamma': gamma,
        'black_level': black_level,
        'black_level_lift': beta,
        'reference_white': ref_white,
        'surround_luminance': surround_luminance
    }

def main():
    """
    Main function that handles command-line arguments and runs the appropriate calculation mode.

    Calculation Modes:
        nominal  : Calculates System Gamma and HDR Reference White for a standard display setup
        extended : Calculates nominal luminance of 100% and the extended 109% luminance value

    Required Parameters:
        --peak      : Display peak luminance in cd/m² (default: 1000)
        --surround  : Surround luminance in cd/m² (default: 5)
        --mode      : Calculation mode (default: nominal)

    Optional Parameters:
        --black   : Enable black level lift with specified value in cd/m² (default: 0.0)
        --pluge   : Calculate PLUGE test signal values according to ITU-R BT.814
        --explain : Show detailed explanations of calculation steps

    Usage examples:
        python3 hlg_calculator.py --peak 1000 --surround 5 --mode nominal --explain
        python3 hlg_calculator.py --peak 1000 --surround 5 --pluge --explain
        python3 hlg_calculator.py --peak 1000 --surround 5 --black 0.005 --pluge --mode nominal --explain
        python3 hlg_calculator.py --peak 1000 --surround 5 --mode extended --explain
    """
    # Command-line argument parser setup
    parser = argparse.ArgumentParser(description="HLG Display Adaptation Calculator")

    parser.add_argument('--peak', type=float, default=1000.0,
                        help="Peak luminance in cd/m² (default: 1000)")
    parser.add_argument('--surround', type=float, default=5.0,
                        help="Surround luminance in cd/m² (default: 5)")
    parser.add_argument('--black', nargs='?', const=0.005, type=float, default=0.0,
                        help="Enable black level lift with specified value in cd/m² (default when flag used: 0.005)")
    parser.add_argument('--pluge', action='store_true',
                        help="Calculate PLUGE test signal values according to ITU-R BT.814")
    parser.add_argument('--mode', choices=['nominal', 'extended'],
                        default='nominal', help="Gamma calculation mode (default: nominal)")
    parser.add_argument('--explain', action='store_true',
                        help="Show detailed explanations of calculation steps")

    args = parser.parse_args()

    # Set explanation mode based on command-line flag
    global EXPLAIN
    EXPLAIN = args.explain

    # Calculation modes output
    if args.mode == 'nominal':
        # Calculate Nominal Range
        result = nominal_range(args.peak, args.surround, args.black)

        print("\nNominal Range Calculation:")
        print(f"  Nominal 100% luminance:    {result['100%_luminance']:.0f} cd/m²")
        print(f"  Surround luminance:        {result['surround_luminance']:.0f} cd/m²")
        if args.black > 0.0:
            print(f"  Black level luminance:     {result['black_level']:.4f} cd/m²")
            print(f"  Black level lift (β):      {result['black_level_lift']:.4f}")
        print(f"  System gamma:              {result['system_gamma']:.2f}")
        print(f"  HDR Reference White (75%): {result['reference_white']:.0f} cd/m²\n")

    elif args.mode == 'extended':
        # Calculate 109% Extended Range luminance based on 100% nominal range luminance
        result = extended_range(args.peak, args.surround, args.black)
        
        print("\nExtended Range Calculation:")
        print(f"  Nominal 100% luminance:    {result['100%_luminance']:.0f} cd/m²")
        print(f"  Extended 109% luminance:   {result['109%_luminance']:.0f} cd/m²")
        print(f"  Surround luminance:        {result['surround_luminance']:.0f} cd/m²")
        if args.black > 0.0:
            print(f"  Black level luminance:     {result['black_level']:.4f} cd/m²")
            print(f"  Black level lift (β):      {result['black_level_lift']:.4f}")
        print(f"  System gamma:              {result['system_gamma']:.2f}")
        print(f"  HDR Reference White (75%): {result['reference_white']:.0f} cd/m²\n")

    # Calculate PLUGE values if the pluge argument is used
    if args.pluge:
        pluge_result = calculate_pluge_values(args.peak, args.surround, args.black)

        # Print PLUGE results
        print("\nPLUGE Test Signal Luminance Values (10-bit/12-bit)")
        print("According to ITU-R BT.814 Table 3:")
        print(f"  Higher level (399/1596):     {pluge_result['higher_level']:.4f} cd/m²")
        print(f"  Slightly lighter (80/320):   {pluge_result['slightly_lighter_level']:.4f} cd/m²")
        print(f"  Black level (64/256):        {pluge_result['black_level_display']:.4f} cd/m²")
        print(f"  Slightly darker (48/192):    {pluge_result['slightly_darker_level']:.4f} cd/m²")
        print()
        print("  Note: Values correspond to the HDR code values in ITU-R BT.814 Table 3")
        print("        To correctly set up a display, the \"slightly darker\" level should be")
        print("        just below visibility threshold while \"slightly lighter\" remains visible.\n")

if __name__ == "__main__":
    main()
