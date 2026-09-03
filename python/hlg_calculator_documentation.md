# HLG Display Adaptation Calculator Documentation

Documentation for `hlg_calculator.py`, with detailed explanations for each function's purpose, mathematical basis and how the different functions and formulas are used.

This HLG Display Adaptation Calculator implements the Hybrid Log-Gamma (HLG) system as specified in ITU-R BT.2100, with calculations for system gamma and HDR Reference White following *ITU-R BT.2100-3 Note 5f*, *ITU-R BT.2390 Section 6.2* and *ITU-R BT.2408 Section 2.1* respectively.

The intention is to help users understand more about how the Hybrid Log-Gamma (HLG) system works.

The idea for this tool was inspired by work related to [EBU Recommendation R 167](https://tech.ebu.ch/publications/r167) — *Reference Monitors: Predefined Modes for HLG*.

## Features and Usage
```
Calculation Modes:
    nominal      : Calculates system gamma and HDR Reference White for a standard display setup
    extended     : Calculates nominal luminance of 100% and the extended 109% luminance value
    extendedpeak : Calculates the nominal peak luminance so the full extended range (0-109%)
                   fits within --peak; lowering the effective nominal peak this way is the
                   ITU-R BT.2100-3 Note 5f contrast-control adjustment

Required Parameters:
    --peak     : Nominal peak luminance in cd/m² (default: 1000);
                 in extendedpeak mode, the Extended Peak Luminance
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
    python3 hlg_calculator.py --peak 2000 --surround 5 --mode extendedpeak --explain
```

## Signal Chain, Modes and Function Relationships

```
1. Calculate System Gamma (γ)
    ├── calculate_system_gamma()
    ├── Uses ITU-R BT.2100-3 Note 5f for peak luminance
    ├── Uses ITU-R BT.2390 for surround luminance

2. Calculate Black Level Lift (if enabled)
    ├── calculate_black_level_lift()
    ├── Uses γ from calculate_system_gamma()

3. Apply OOTF (Perceptual Luminance Adjustment)
    ├── hlg_reference_ootf()
    ├── Uses γ from calculate_system_gamma()

4. Convert Signal to Display Light (HLG EOTF)
    ├── hlg_reference_eotf()
    ├── Uses β from calculate_black_level_lift()
    ├── Maps HLG signal (E') → Display luminance (Fd)

5. Calculate HDR Reference White (75% HLG)
    ├── calculate_reference_white()
    ├── Uses hlg_reference_eotf() to calculate 75% white luminance

6. Calculate PLUGE values (if enabled)
    ├── calculate_pluge_values()
    ├── Uses hlg_reference_eotf() to check black level visibility
```

```shell
nominal_range()
    ├── calculate_system_gamma()       # Calculate γ
    ├── calculate_black_level_lift()   # Calculate β (if enabled)
    ├── calculate_reference_white()    # Calculate 75% Reference White
    └── Print Output
```

```shell
extended_range()
      ├── calculate_system_gamma()       # Calculate γ
      ├── calculate_black_level_lift()   # Calculate β (if enabled)
      ├── hlg_reference_eotf()           # Calculate luminance at 109% signal
      ├── calculate_reference_white()    # Calculate 75% Reference White
      └── Print Output
```

```shell
extended_peak_range()
      ├── compute_peak_109()               # 109% luminance for a trial Lw
      │     ├── calculate_system_gamma()   # γ for the trial Lw
      │     ├── calculate_black_level_lift()  # β (if enabled)
      │     └── hlg_reference_eotf()       # EOTF output at the 109% signal level
      ├── largest_fitting_lw()             # Bisection; second pass from Lw = 400
      └── extended_range()                 # Final parameters for the solved Lw
```

## Order of operations in this script

The HLG Display Adaptation Calculator implements a comprehensive system where functions build upon each other in a logical progression:

1. **System Parameters → Transform Functions**
   - System gamma and black level lift are calculated first
   - These parameters customize the HLG transforms to the specific display and environment
   - The transform functions (EOTF, OOTF) then use these parameters to render content correctly

2. **HLG Signal Chain Integration**
   - The EOTF integrates multiple steps:
     1. Black level lift application
     2. Inverse OETF to convert signal to scene light
     3. OOTF to apply gamma and convert to display light
   - This matches the signal chain specified in ITU-R BT.2100-3

3. **Black Level Lift and PLUGE**
   - Black level lift ($\beta$) and PLUGE are closely integrated:
     - PLUGE test signals allow proper setting of display black level ($L_B$)
     - The measured black level is used to calculate the black level lift parameter
     - The lift parameter ensures proper rendering of dark details
   - ITU-R BT.2100-3 Note 5h specifically states: "Values below 0.0 should not be clipped in reference displays (even though they represent 'negative' light) to allow the black level of the signal ($L_B$) to be properly set using test signals known as 'PLUGE'."

4. **Calculation Modes**

    - The three calculation modes use the functions to calculate the correct system gamma, Reference White, black level.

      - **nominal**      : Calculates system gamma and HDR Reference White for a standard display setup
      - **extended**     : Calculates nominal luminance of 100% and the extended 109% luminance value
      - **extendedpeak** : Calculates the nominal peak luminance so the full extended range (0–109%) fits within a given Extended Peak Luminance


## 1. HLG OETF and Inverse OETF

This section implements the Hybrid Log-Gamma Reference Opto-Electronic Transfer Function (OETF) and its inverse as specified in ITU-R BT.2100-3 Table 5. These functions provide the fundamental transforms between scene linear light and non-linear HLG signal values.

### HLG Reference OETF

```python
def hlg_reference_oetf(E):
    """
    HLG Opto-Electronic Transfer Function (OETF) as defined in ITU-R BT.2100-3 Table 5.
    
    Transforms scene linear light to the non-linear HLG signal.
    
    Args:
        E (float): Scene linear light normalized to the range [0:1]
        
    Returns:
        float: Non-linear HLG signal value E' in the range [0:1]

    Where diffuse white and an 18% grey card land on this curve is set by the
    camera exposure ahead of the OETF (BT.2100-3 Note 5b), not by the OETF.
        
    Reference: ITU-R BT.2100-3, Table 5, HLG Reference OETF
    """
    # Implement the piecewise HLG OETF function
    if E <= 1/12:
        return math.sqrt(3 * E)
    else:
        return HLG_a * math.log(12 * E - HLG_b) + HLG_c
```

The OETF is the camera encoding function that converts scene linear light to non-linear HLG signal values. It uses a piecewise function with a square-root curve for dark regions (E ≤ 1/12) and a logarithmic curve for bright regions (E > 1/12). The constants HLG_a, HLG_b, and HLG_c are precisely defined in the ITU-R BT.2100-3 standard. Where diffuse white and an 18% grey card land on the curve is an exposure decision made ahead of the OETF (BT.2100-3 Note 5b).

### HLG Inverse OETF

```python
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
        float: Scene linear light (E), normalized so that E = 1 corresponds
               to the 100% signal; exceeds 1.0 for super-white input (E' > 1).
               A value of 1/12 corresponds to the transition point (E'=0.5).
    
    Reference: ITU-R BT.2100-3, Table 5, Hybrid Log-Gamma (HLG) system,
               "inverse OETF" component of the "Reference HLG EOTF" definition.
    """
    # Ensure signal is non-negative per BT.2100-3 Table 5
    E_prime = max(0.0, E_prime)
    
    # Apply the piecewise inverse OETF
    if E_prime <= 0.5:
        # Square law segment (E ≤ 1/12 in original OETF)
        E = (E_prime ** 2) / 3.0
    else:
        # Logarithmic segment (E > 1/12 in original OETF)
        # Calculate intermediate values
        exp_term = math.exp((E_prime - HLG_c) / HLG_a)
        E = (exp_term + HLG_b) / 12.0
    
    return E
```

The Inverse OETF converts non-linear HLG signals back to scene linear light. It's the mathematical inverse of the OETF and is a critical component of the EOTF (display rendering) process. The 0.5 threshold corresponds to the 1/12 threshold in the forward OETF.

## 2. HLG System Gamma and Black Level Lift

This section calculates the critical display adaptation parameters specified in ITU-R BT.2100-3. The system gamma calculation implements the formula in Note 5f, extended with the surround luminance adjustment from ITU-R BT.2390 Section 6.2. The black level lift parameter ($\beta$) is calculated according to the formula specified in Table 5 of ITU-R BT.2100-3.

### Calculate System Gamma

```python
def calculate_system_gamma(Lw, Ls, step_number=None):
    """
    Calculate the system gamma (γ) according to formulas specified in:
    - ITU-R BT.2100-3 Note 5f
    - ITU-R BT.2390 Section 6.2
    """
    # Constants
    gamma_ref = 1.2  # γ_ref
    kappa = 1.111    # κ
    mu = 0.98        # μ
    L_ref = 1000     # L_ref
    Ls_ref = 5       # Ls_ref
    
    # Select the appropriate formula based on nominal peak luminance
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
    mu_factor = mu ** math.log2(Ls / Ls_ref)
    
    # Apply surround adjustment to system gamma
    gamma = base_gamma * mu_factor
    
    return gamma, {
        'gamma_ref': gamma_ref,
        'kappa': kappa,
        'mu': mu,
        'L_ref': L_ref,
        'Ls_ref': Ls_ref,
        'base_gamma': base_gamma,
        'kappa_factor': kappa_factor,
        'mu_factor': mu_factor
    }
```

This function calculates the system gamma ($\gamma$) based on two adaptation conditions:
1. The display's nominal peak luminance ($L_W$)
2. The viewing environment's surround luminance ($L_{surround}$)

System gamma is a critical parameter that adapts the HLG display response to maintain perceptual consistency across different display capabilities and viewing conditions:
- Higher nominal peak luminance requires increased gamma ($\kappa > 1$)
- Brighter viewing environments require decreased gamma ($\mu < 1$)

#### System Gamma Formulas

Two different formulas are specified in ITU-R BT.2100-3 Note 5f, each optimized for a different nominal peak luminance range:

- For $400 \text{ cd/m}^2 \leq L_W \leq 2000 \text{ cd/m}^2$: 
  $$\gamma = 1.2 + 0.42 \log_{10}\left(\frac{L_W}{1000}\right)$$
  
- For $L_W < 400 \text{ cd/m}^2$ or $L_W > 2000 \text{ cd/m}^2$: 
  $$\gamma = 1.2 \times \kappa^{\log_2\left(\frac{L_W}{1000}\right)} \text{ where } \kappa = 1.111$$

According to ITU-R BT.2100-3 Note 5f:
> "For displays with nominal peak luminance ($L_W$) other than 1 000 cd/m², or where the effective nominal peak luminance is adjusted through the use of a contrast control, the system gamma value should be adjusted according to the formula below and may be rounded to three significant digits:"
>
> $\gamma = 1.2 + 0.42 \log_{10}(L_W/1000)$
>
> "For applications outside of the usual production monitoring range of $L_W$ equal to 400 cd/m² to 2 000 cd/m², the following extended range formula should be used:"
>
> $\gamma = 1.2 \times \kappa^{\log_2(L_W/1000)}$ where $\kappa = 1.111$

Both formulas are centered around the reference point of 1000 cd/m², where they produce identical results with a system gamma of $\gamma = 1.2$.

#### Surround Luminance Adjustment

ITU-R BT.2390 Section 6.2 describes how system gamma should be reduced when the viewing environment is brighter than the 5 cd/m² reference surround. The Report presents two models for this; the calculator implements the multiplicative model with the surround luminance compensation factor ($\mu$). In both cases above, this additional adjustment for surround luminance should be applied:

$$\gamma = \gamma \times \mu^{\log_2\left(\frac{L_{surround}}{L_{surround\ ref}}\right)}$$

where $\mu = 0.98$ and $L_{surround\ ref} = 5 \text{ cd/m}^2$

#### Complete Formula

The complete extended formula from ITU-R BT.2390 Section 6.2 is:

$$\gamma = \gamma_{ref} \times \kappa^{\log_2\left(\frac{L_W}{L_{ref}}\right)} \times \mu^{\log_2\left(\frac{L_{surround}}{L_{surround\ ref}}\right)}$$

Where:
- $\gamma_{ref} = 1.2$ (reference gamma at 1000 cd/m² in reference environment)
- $\kappa = 1.111$ (peak luminance adjustment factor)
- $\mu = 0.98$ (surround luminance adjustment factor)
- $L_{ref} = 1000 \text{ cd/m}^2$ (reference nominal peak luminance)
- $L_{surround\ ref} = 5 \text{ cd/m}^2$ (reference surround luminance)

#### When to use each formula?

| Condition | Recommended Formula |
|-----------|---------------------|
| $L_W$ changes, $L_{surround}$ fixed (e.g., reference environment with $L_{surround} = 5 \text{ cd/m}^2$) | ITU-R BT.2100-3 Note 5f |
| Both $L_W$ and $L_{surround}$ change (e.g., viewing in varying ambient conditions) | ITU-R BT.2390 Section 6.2 |
| $L_W$ between 400-2000 cd/m², $L_{surround}$ changes | Basic formula with surround adjustment |
| $L_W$ outside 400-2000 cd/m², $L_{surround}$ changes | Extended formula with surround adjustment |

#### Real-world Applications

In real-world applications, both the nominal peak luminance and viewing conditions vary. The human visual system adapts to both the brightness of the display ($L_W$) and the ambient lighting ($L_{surround}$). The formulas specified in ITU-R BT.2100-3 Note 5f and ITU-R BT.2390 Section 6.2 ensures perceptual consistency across different displays and viewing environments.

### Calculate Black Level Lift

```python
def calculate_black_level_lift(L_B, Lw, gamma, step_number=None):
    """
    Calculate the black level lift (β) according to ITU-R BT.2100-3.
    
    The black level lift parameter adjusts signal mapping to accommodate non-zero
    black levels in real displays, ensuring consistent rendering across displays
    with different black level capabilities.
    
    Per ITU-R BT.2100-3 Table 5, black level lift is defined as:
        β = √(3(L_B/Lw)^(1/γ))

    Where:
    - L_B is the display luminance for black in cd/m²
    - Lw is nominal peak luminance of the display in cd/m²
    - γ is the system gamma
    
    This parameter ensures proper rendering of the signal's black level on displays
    where absolute zero luminance isn't achievable. It modifies the signal mapping
    to maintain perceptual consistency across different display technologies.
    
    Args:
        L_B (float): Display luminance for black in cd/m².
        Lw (float): Nominal peak luminance of the display in cd/m².
        gamma (float): System gamma calculated for the display and viewing environment.
        step_number (int, optional): Step number for explanation purposes.
        
    Returns:
        float: Black level lift factor (β).
    
    Reference: ITU-R BT.2100-3, Table 5, HLG Reference EOTF, definition of β.
    """
    # Early return if black level is 0 (ideal display)
    if L_B <= 0:
        return 0.0
        
    # Calculate black level lift: β = √(3(L_B/Lw)^(1/γ))
    beta = math.sqrt(3 * ((L_B / Lw) ** (1 / gamma)))
    
    return beta
```

Black level lift ($\beta$) compensates for non-zero black levels in real displays. It adjusts signal mapping to ensure consistent rendering across displays with different black level capabilities.

The formula $\beta = \sqrt{3(L_B/L_W)^{1/\gamma}}$ ensures proper rendering of the signal's black level on displays where absolute zero luminance isn't achievable.

The black level lift factor as specified in ITU-R BT.2100-3 is calculated as:

$$\beta = \sqrt{3 \times \left(\frac{L_B}{L_W}\right)^{\frac{1}{\gamma}}}$$

Where:

- $L_B$ is the display luminance for black in cd/m²
- $L_W$ is nominal peak luminance of the display in cd/m²
- $\gamma$ is the system gamma calculated for the specific display and viewing environment

## 3. HLG Reference OOTF

This section implements the HLG Reference Opto-Optical Transfer Function (OOTF) as defined in ITU-R BT.2100-3 Table 5. The OOTF maps scene linear RGB to display linear RGB, applying System Gamma to the luminance component only as specified in the standard. This approach preserves the chromaticity of the scene as imaged by the camera.


```python
def hlg_reference_ootf(r_s, g_s, b_s, gamma, alpha=1.0, explain_label=None):
    """
    HLG Reference Opto-Optical Transfer Function (OOTF) as defined in ITU-R BT.2100-3 Table 5.
    
    Maps scene linear RGB to display linear RGB with appropriate gamma adjustment.
    The HLG OOTF applies the system gamma to the luminance component only,
    preserving consistent color appearance through the rendering chain.
    
    As specified in ITU-R BT.2390 Section 6.2:
    "Instead of the current SDR practice of applying a gamma curve independently
    to each colour component, for HDR it should be applied to the luminance alone."
    
    The HLG OOTF calculates a single luminance adjustment factor (Y_S^(γ-1))
    and applies it uniformly to all three color components. This preserves the
    ratios between R, G, and B, maintaining consistent color appearance across
    different luminance levels.
    
    Args:
        r_s, g_s, b_s (float): Scene linear RGB components [0:1]
        gamma (float): System gamma calculated based on display and environment
        alpha (float, optional): User gain α in cd/m², which BT.2100-3 Table 5
                                 defines as representing Lw, the nominal peak
                                 luminance; the default 1.0 returns display
                                 light normalized to the nominal peak
        explain_label (str, optional): Label for explanation output
    
    Returns:
        tuple: Display linear RGB components (r_d, g_d, b_d), in cd/m² when
               alpha is Lw
    
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
    
    return r_d, g_d, b_d
```

The OOTF maps scene linear RGB to display linear RGB with appropriate gamma adjustment. It implements a critical conceptual difference for HDR:

"Instead of the current SDR practice of applying a gamma curve independently to each colour component, for HDR it should be applied to the luminance alone." - ITU-R BT.2390 Section 6.2

The OOTF calculates a single luminance adjustment factor ($Y_S^{(\gamma-1)}$) and applies it uniformly to all three color components, preserving the ratios between R, G, and B and maintaining consistent color appearance.

## 4. HLG Reference EOTF (using Inverse OETF + OOTF)

This section implements the complete HLG Reference Electro-Optical Transfer Function (EOTF) as specified in ITU-R BT.2100-3 Table 5. The implementation follows the standard's definition of the EOTF as a three-step process: applying black level lift, converting to scene light via the inverse OETF, and applying the OOTF to produce display light.

```python
def hlg_reference_eotf(r_prime, g_prime, b_prime, gamma, Lw, L_B=0.0, beta=None, explain_label=None):
    """
    Implements the full HLG Electro-Optical Transfer Function (EOTF) with OOTF 
    for RGB signals as specified in ITU-R BT.2100-3 Table 5, with optional
    black level lift support.
    
    The HLG Reference EOTF consists of these conceptual steps:
    1. Apply black level lift if required: E' → max(0,(1-β)E'+β)
    2. Inverse OETF: Convert non-linear R', G', B' to scene linear R, G, B
    3. Apply OOTF: F_D = α·Y_S^(γ-1)·E, where:
       - Y_S is the scene luminance: 0.2627R + 0.6780G + 0.0593B
       - γ is the variable system gamma
       - α is the user gain in cd/m², representing Lw, so F_D is the
         display luminance in cd/m² (BT.2100-3 Table 5)
    
    The complete chain is expressed as: F_D = OOTF[OETF⁻¹[max(0,(1-β)E'+β)]]
    
    Args:
        r_prime, g_prime, b_prime (float): Non-linear HLG signal values [0:1+]
        gamma (float): System gamma calculated according to ITU-R BT.2100-3 Note 5f
                       Typically in range 1.0-1.5, based on nominal peak luminance
        Lw (float): Nominal peak luminance of the display in cd/m²; the
                    OOTF's user gain α and the reference for black level lift
        L_B (float, optional): Display luminance for black in cd/m²
        beta (float, optional): Precalculated black level lift factor to reuse
                               If None, will be calculated from L_B and Lw
        explain_label (str, optional): Label for explanation output
    
    Returns:
        tuple: Display luminance (r_d, g_d, b_d) in cd/m²
    
    Reference: ITU-R BT.2100-3, Table 5, "HLG Reference EOTF" definition.
    """
    # Step 1: Apply black level lift if required (L_B > 0)
    black_level_enabled = L_B > 0.0 or (beta is not None and beta > 0)
    
    if black_level_enabled:
        # Use provided beta or calculate it
        if beta is None:
            beta = calculate_black_level_lift(L_B, Lw, gamma)
        
        # Apply black level lift to each component per BT.2100-3 formula:
        # F_D = EOTF[max(0,(1-β)E'+β)]
        r_lifted = max(0, (1 - beta) * r_prime + beta)
        g_lifted = max(0, (1 - beta) * g_prime + beta)
        b_lifted = max(0, (1 - beta) * b_prime + beta)
    else:
        # No black level lift
        r_lifted = r_prime
        g_lifted = g_prime
        b_lifted = b_prime
    
    # Step 2: Convert non-linear signal to scene linear light using inverse OETF
    # Only pass the explain_label for achromatic signals to avoid redundant explanations
    explain_inverse_oetf = explain_label and r_lifted == g_lifted == b_lifted
    
    r_s = hlg_inverse_oetf(r_lifted, explain_label if explain_inverse_oetf else None)
    g_s = hlg_inverse_oetf(g_lifted)
    b_s = hlg_inverse_oetf(b_lifted)
    
    # Step 3: Apply the OOTF with α = Lw to get display luminance in cd/m²
    r_d, g_d, b_d = hlg_reference_ootf(r_s, g_s, b_s, gamma, alpha=Lw, explain_label=explain_label)
    
    return r_d, g_d, b_d
```

The EOTF is the complete display rendering function that maps non-linear HLG signals to display light. It combines three conceptual steps:
1. Apply black level lift if required: $E' → \max(0,(1-\beta)E'+\beta)$
2. Convert non-linear signals to scene linear using the inverse OETF
3. Apply the OOTF with α = L<sub>W</sub> to get display luminance in cd/m²

This function is the central integration point of the HLG system, bringing together black level compensation, inverse signal mapping, and the critical rendering transforms.

### HLG EOTF — End-to-End Signal Chain
The complete HLG EOTF with black level lift can be conceptually represented as:

$$F_D = \mathrm{OOTF}[\mathrm{OETF}^{-1}[\max(0,(1-\beta)E'+\beta)]]$$

Where:

- $F_D$ is the display-referred linear light
- $\text{OETF}^{-1}$ is the inverse OETF (Opto-Electronic Transfer Function)
- $\text{OOTF}$ is the Opto-Optical Transfer Function
- $E'$ is the non-linear HLG signal value
- $\beta$ is the black level lift factor

The OOTF applies gamma to the luminance component only, which can be represented as:

$$F_D = \alpha \cdot Y_S^{(\gamma-1)} \cdot E$$

Where:

- $Y_S$ is the scene luminance calculated using BT.2020 coefficients:

$$Y_S = 0.2627R_S + 0.6780G_S + 0.0593B_S$$

- $\gamma$ is the system gamma calculated using the complete formula
- $\alpha$ is the user gain in cd/m², which *ITU-R BT.2100-3 Table 5* defines as representing $L_W$, the nominal peak luminance, so $F_D$ is display luminance in cd/m²
- $E$ represents the scene-referred linear RGB components

For achromatic signals (where R=G=B), this simplifies to:

$$F_D = \alpha \cdot E^{\gamma} = L_W \cdot E^{\gamma}$$

## 5. HDR Reference White and PLUGE

This section calculates the HDR Reference White signal level as defined in ITU-R BT.2408 Section 2.1 and the PLUGE test signal values for proper black level calibration according to ITU-R BT.814.

### Calculate Reference White

```python
def calculate_reference_white(Lw, gamma, L_B=0.0, beta=None, step_number=None):
    """
    Calculate the luminance of HDR Reference White (75% signal level).

    Based on ITU-R BT.2408 Section 2.1, HDR Reference White is defined as:
    "the nominal signal level obtained from an HDR camera and a 100% reflectance 
    white card resulting in a nominal luminance of 203 cd/m² on a PQ display or 
    on an HLG display that has a nominal peak luminance capability of 1,000 cd/m²."
    
    For the 1000 cd/m² reference display with γ = 1.2 this function returns
    203.15 cd/m², which rounds to the 203 cd/m² of ITU-R BT.2408 Table 1. The
    75% signal is BT.2408's rounded nominal; 203 cd/m² exactly is 74.99%.
    
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
    
    # The EOTF returns display luminance in cd/m² (α = Lw), with black level
    # lift applied if enabled
    ref_white, _, _ = hlg_reference_eotf(signal_75, signal_75, signal_75, gamma,
                                         Lw=Lw, L_B=L_B, beta=beta, explain_label="75%")
    
    return ref_white
```

This function calculates the luminance of HDR Reference White (75% signal level) based on ITU-R BT.2408 Section 2.1:

"HDR Reference White is the nominal signal level obtained from an HDR camera and a 100% reflectance white card resulting in a nominal luminance of 203 cd/m² on a PQ display or on an HLG display that has a nominal peak luminance capability of 1,000 cd/m²."

For the 1000 cd/m² reference display with γ = 1.2 this function returns 203.15 cd/m², which rounds to the 203 cd/m² of ITU-R BT.2408 Table 1. The 75% signal is BT.2408's rounded nominal; 203 cd/m² exactly is 74.99%.

### Calculate PLUGE Values

```python
def calculate_pluge_values(peak_luminance, surround_luminance, black_level=0.0):
    """
    Calculate PLUGE test signal levels for HLG displays according to ITU-R BT.814.
    
    This function calculates the luminance values corresponding to the 10-bit and 12-bit
    code values specified in ITU-R BT.814 Table 3 for HDR displays. It uses the HLG EOTF
    with system gamma adapted to the viewing environment.
    
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
            - 'system_gamma': The calculated system gamma
            - 'black_level_lift': The calculated black level lift factor (β)
            - 'input_black_level': Original black level input parameter
    """
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
    
    # 12-bit values (used for verification output in explain mode)
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
    
    # Calculate normalized values from 10-bit codes
    higher_level_normalized = (higher_level_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    black_level_normalized = (black_level_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    slightly_lighter_normalized = (slightly_lighter_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    slightly_darker_normalized = (slightly_darker_code_10bit - black_10bit) / (peak_10bit - black_10bit)
    
    # Calculate luminance values using the HLG EOTF for each PLUGE level
    # (α = Lw, so the EOTF returns cd/m² directly)
    higher_level_luminance, _, _ = hlg_reference_eotf(higher_level_normalized, higher_level_normalized, higher_level_normalized, 
        gamma, Lw=peak_luminance, L_B=black_level, beta=beta)
    
    black_level_display_luminance, _, _ = hlg_reference_eotf(black_level_normalized, black_level_normalized, black_level_normalized, 
        gamma, Lw=peak_luminance, L_B=black_level, beta=beta)
    
    slightly_lighter_luminance, _, _ = hlg_reference_eotf(slightly_lighter_normalized, slightly_lighter_normalized, slightly_lighter_normalized, 
        gamma, Lw=peak_luminance, L_B=black_level, beta=beta)
    
    slightly_darker_luminance, _, _ = hlg_reference_eotf(slightly_darker_normalized, slightly_darker_normalized, slightly_darker_normalized, 
        gamma, Lw=peak_luminance, L_B=black_level, beta=beta)
    
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
```

This function calculates PLUGE test signal levels according to ITU-R BT.814 for HLG displays. It calculates the luminance values for key reference points:
- Higher level (code value 399/1596)
- Black level (code value 64/256)
- Slightly lighter level (code value 80/320)
- Slightly darker level (code value 48/192)

These values are essential for proper black level adjustment of displays. The slightly darker level should be just below visibility threshold while slightly lighter level should remain visible.

## 6. Mode-specific Calculation Functions

### Nominal Range

```python
def nominal_range(peak_luminance, surround_luminance, black_level=0.0):
    """
    Calculate standard parameters for an HLG display.
    
    This function calculates the essential parameters for an HLG display setup,
    including system gamma adaptation based on viewing environment, optional
    black level lift, and HDR Reference White level.

    Args:
        peak_luminance (float): Nominal peak luminance of the display in cd/m².
        surround_luminance (float): Luminance of the display surround in cd/m².
        black_level (float, optional): Display luminance for black in cd/m².

    Returns:
        dict: A dictionary containing:
            - '100%_luminance' (float): Nominal 100% luminance.
            - 'system_gamma' (float): Calculated system gamma.
            - 'black_level' (float): Display luminance for black.
            - 'black_level_lift' (float): Calculated black level lift factor (β).
            - 'reference_white' (float): HDR Reference White in cd/m².
            - 'surround_luminance' (float): Input surround luminance (for reference).
    """
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
    
    # Step 3: Calculate HDR Reference White
    ref_white = calculate_reference_white(peak_luminance, gamma, 
                                        L_B=black_level, beta=beta, 
                                        step_number=ref_white_step)
    
    return {
        '100%_luminance': peak_luminance,
        'system_gamma': gamma,
        'black_level': black_level,
        'black_level_lift': beta,
        'reference_white': ref_white,
        'surround_luminance': surround_luminance
    }
```

This function calculates the standard parameters for an HLG display setup, integrating system gamma adaptation, optional black level lift, and HDR Reference White level.

### Extended Range

```python
def extended_range(peak_luminance, surround_luminance, black_level=0.0,
                   nominal_is_calculated=False, first_step=1):
    """
    Calculate extended range parameters for an HLG display, showing the brightness
    at 109% signal level (super-white) for a given nominal 100% luminance.

    This function provides detailed calculation steps to show how HLG handles signals
    beyond the nominal range, with optional support for black level lift.

    Args:
        peak_luminance (float): Nominal peak luminance of the display in cd/m²
        surround_luminance (float): Luminance of the display surround in cd/m²
        black_level (float, optional): Display luminance for black in cd/m²
        nominal_is_calculated (bool, optional): When True, the explanation
            header labels the nominal peak luminance as calculated rather
            than as an input (used by extended_peak_range(), where Lw is
            solved, not given). The calculations are unaffected.
        first_step (int, optional): Number of the first explanation step.
            extended_peak_range() passes 2 so its step numbering continues
            past the solver's Step 1 instead of restarting.

    Returns:
        dict: A dictionary containing:
            - '100%_luminance' (float): Nominal 100% luminance.
            - '109%_luminance' (float): Calculated extended luminance at 109% signal.
            - 'system_gamma' (float): Calculated system gamma.
            - 'black_level' (float): Display luminance for black.
            - 'black_level_lift' (float): Calculated black level lift factor (β).
            - 'reference_white' (float): HDR Reference White in cd/m².
    """
    # Step 1: Calculate System Gamma
    gamma, _ = calculate_system_gamma(peak_luminance, surround_luminance, step_number=first_step)

    # Step 2: Calculate Black Level Lift (if enabled)
    beta = 0.0
    step_offset = 0
    if black_level > 0.0:
        beta = calculate_black_level_lift(black_level, peak_luminance, gamma, step_number=first_step + 1)
        step_offset = 1

    # Step 3: Calculate 109% signal luminance
    E = SIGNAL_LEVEL_109

    # The EOTF returns display luminance in cd/m² (α = Lw), with black level
    # lift applied if enabled
    luminance_109, _, _ = hlg_reference_eotf(E, E, E, gamma, Lw=peak_luminance,
                                             L_B=black_level, beta=beta,
                                             explain_label="109%")
    
    # Step 4: Calculate HDR Reference White
    ref_white = calculate_reference_white(peak_luminance, gamma, L_B=black_level,
                                        beta=beta, step_number=first_step + 2 + step_offset)
    
    return {
        '100%_luminance': peak_luminance,
        '109%_luminance': luminance_109,
        'system_gamma': gamma,
        'black_level': black_level,
        'black_level_lift': beta,
        'reference_white': ref_white,
        'surround_luminance': surround_luminance
    }
```

This function expands on the nominal range calculations by showing how the HLG system handles signals beyond the nominal range, specifically at the standardized 109% extended-range signal level (super-white).

### Extended Peak

```python
def extended_peak_range(extended_peak, surround_luminance, black_level=0.0):
    """
    Solve for the nominal peak luminance (Lw) whose full extended signal range
    (0-109%) fits within a given Extended Peak Luminance, then calculate the
    extended range parameters for that nominal peak.

    The Extended Peak Luminance is the luminance the 109% super-white signal
    level is mapped to. Lowering the effective nominal peak so that the
    extended range fits within the display capability is the "contrast
    control" adjustment described in ITU-R BT.2100-3 Note 5f.

    Args:
        extended_peak (float): Extended Peak Luminance in cd/m² — the
                               luminance the 109% signal level is mapped to.
        surround_luminance (float): Luminance of the display surround in cd/m².
        black_level (float, optional): Display luminance for black in cd/m².

    Returns:
        dict: The extended_range() result dictionary for the solved nominal
              peak luminance. '109%_luminance' is the actual 109% luminance
              of the solved nominal peak; within the no-solution window it
              is slightly below the requested Extended Peak Luminance.
    """
    if extended_peak <= 0:
        raise ValueError("Extended Peak Luminance must be positive")

    def compute_peak_109(lw):
        gamma, _ = calculate_system_gamma(lw, surround_luminance)
        beta = 0.0
        if black_level > 0.0:
            beta = calculate_black_level_lift(black_level, lw, gamma)
        E = SIGNAL_LEVEL_109
        peak_109, _, _ = hlg_reference_eotf(E, E, E, gamma, Lw=lw,
                                            L_B=black_level, beta=beta)
        return peak_109

    # One shared iteration limit across both bisection passes
    iterations = 0

    def largest_fitting_lw(lo, hi):
        nonlocal iterations
        while (hi - lo) / hi > 1e-7 and iterations < 100:
            iterations += 1
            mid = (lo + hi) / 2
            if compute_peak_109(mid) <= extended_peak:
                lo = mid
            else:
                hi = mid
        return lo

    # EOTF[1.0902] ratio is > 1 and < 3 for any reachable gamma, so the
    # solution always lies between extended_peak / 3 and extended_peak.
    lw = largest_fitting_lw(extended_peak / 3, extended_peak)

    # γ steps down across Lw = 400 (Note 5f formula switch), so a fitting Lw
    # just below 400 can coexist with a larger one at or above it; prefer
    # the larger so mode swaps round-trip.
    if lw < 400 and compute_peak_109(400) <= extended_peak:
        lw = largest_fitting_lw(400, extended_peak)

    return extended_range(lw, surround_luminance, black_level,
                          nominal_is_calculated=True, first_step=2)
```

This function inverts the extended range calculation: instead of taking the nominal peak luminance as an input, it solves for the nominal $L_W$ whose full extended signal range (0–109%) fits within a given Extended Peak Luminance. Lowering the effective nominal peak this way is the contrast-control adjustment described in *ITU-R BT.2100-3 Note 5f* ("where the effective nominal peak luminance is adjusted through the use of a contrast control").

An iterative solver is required because system gamma depends on the nominal peak luminance (*ITU-R BT.2100-3 Note 5f*), the EOTF output at the 109% signal level depends on that gamma, and the 109% luminance is the product of the nominal peak and the EOTF output — so the relationship cannot be inverted in closed form.

Bisection is used rather than fixed-point iteration because Note 5f's switch between its two gamma formulas makes $\gamma(L_W)$ discontinuous at 400 and 2000 cd/m² (the two formulas are "virtually identical" only within the 400–2000 cd/m² monitoring range, per *ITU-R BT.2390-12 Section 6.2*). The discontinuities create a small window of Extended Peak Luminance values (≈3855–3869 cd/m² at reference surround and zero black level) with no exact solution and another (≈667–671 cd/m²) with two. Bisection relies only on the 109% luminance increasing with $L_W$ within each formula branch and converges to the largest fitting $L_W$; a second pass from the $L_W$ = 400 boundary covers the two-solution window, preferring the larger solution so that a nominal peak configuration round-trips through this mode.

## References

### ITU Recommendations & Reports

- [ITU-R BT.2100](https://www.itu.int/rec/R-REC-BT.2100)\
  Image parameter values for high dynamic range television for use in production and international programme exchange

- [ITU-R BT.2390](https://www.itu.int/pub/R-REP-BT.2390)\
  High dynamic range television for production and international programme exchange

- [ITU-R BT.2408](https://www.itu.int/pub/R-REP-BT.2408)\
  Guidelines for operational practices in high dynamic range television production

- [ITU-R BT.814](https://www.itu.int/rec/R-REC-BT.814)\
  Specifications of PLUGE test signals and alignment procedures for setting of brightness and contrast of displays

### EBU Recommendations

- [EBU R 103](https://tech.ebu.ch/publications/r103)\
  Video Signal Tolerance in Digital Television Systems

- [EBU R 167](https://tech.ebu.ch/publications/r167)\
  Reference Monitors: Predefined Modes for HLG

### BBC Research & Development Research Papers:

- [High Dynamic Range Television and Hybrid Log-Gamma](https://www.bbc.co.uk/rd/projects/high-dynamic-range)

- [WHP 283 (2014)](https://www.bbc.co.uk/rd/publications/whitepaper283)\
  Non-linear Opto-Electrical Transfer Functions for High Dynamic Range Television (Tim Borer)

- [WHP 309 (2015)](https://www.bbc.co.uk/rd/publications/whitepaper309)\
  A "Display Independent" High Dynamic Range Television System (T. Borer and A. Cotton)

- [WHP 341 (2019)](https://www.bbc.co.uk/rd/publications/whitepaper341)\
  Modelling Brightness Perception for High Dynamic Range Television (K. C. Noland, M. Pindoria and A. Cotton)

- [WHP 342 (2019)](https://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP342.pdf)\
  A Brightness Measure for High Dynamic Range Television (K. C. Noland and M. Pindoria)

- [WHP 369 (2019)](https://www.bbc.co.uk/rd/publications/display-high-dynamic-range-images-varying-viewing-conditions)\
  Display of High Dynamic Range Images Under Varying Viewing Conditions (Tim Borer)