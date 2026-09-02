/*
Copyright 2026 Thomas Berglund (NRK)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/**
 * HLG Display Adaptation Calculator
 *
 * - HLGCalculator: Core calculation logic
 * - UIController: Manages UI elements and updates
 * - BlackLevelController: Manages black level functionality
 * - FormulaRenderer: Handles KaTeX formula rendering
 * - PresetController: Manages preset buttons
 * - SliderResetController: The resets beside the display sliders
 * - ExtendedRangeController: Handles extended range functionality
 * - HLGAdaptationCalculator: Main application class
 *
 * The graph panel's classes are in hlg_graphs.js, loaded before this
 * file; HLGAdaptationCalculator constructs its GraphPanelController.
 */

//=============================================================================
// CORE CALCULATIONS
//=============================================================================

class HLGCalculator {
    constructor() {
        // Constants from ITU-R BT.2100-3 and BT.2390
        this.GAMMA_REF = 1.2;   // Reference gamma at 1000 cd/m²
        this.KAPPA = 1.111;     // Peak luminance adjustment factor
        this.MU = 0.98;         // Surround luminance adjustment factor
        this.L_REF = 1000;      // Reference nominal peak luminance in cd/m²
        this.LS_REF = 5;        // Reference surround luminance in cd/m²

        // HLG constants from ITU-R BT.2100-3 Table 5
        this.a = 0.17883277;
        this.b = 1 - 4 * this.a;  // Equals 0.28466892
        this.c = 0.5 - this.a * Math.log(4 * this.a);  // Equals 0.55991073

        // BT.2020 luminance coefficients as defined in ITU-R BT.2100-3 Table 5
        this.Y_R = 0.2627;  // Red luminance coefficient
        this.Y_G = 0.6780;  // Green luminance coefficient 
        this.Y_B = 0.0593;  // Blue luminance coefficient

        // 109% extended range signal level per EBU R 103 (1.090182648401826),
        // derived from the ITU-R BT.2100-3 Table 9 10-bit narrow-range code values
        this.SIGNAL_LEVEL_109 = (1019 - 64) / (940 - 64);
    }

    /**
     * HLG Reference OETF: Converts scene linear light (E, normalized
     * to [0:1] with E = 1 at the 100% signal) to the non-linear
     * signal value E'.
     *
     * Table 5 defines the OETF over 0 ≤ E ≤ 1, the lower bound being
     * part of the first branch condition, so a negative E is outside
     * the domain and the Recommendation prescribes nothing for it.
     * Clamping it to 0 keeps √(3E) from returning NaN and breaking
     * the plotted curve; Table 5's own max(0, ...) sits in the EOTF,
     * on the lifted signal, and is not part of this OETF.
     *
     * The upper bound is deliberately not enforced. Callers ask for
     * the super-whites: the OETF view's dashed segment runs the log
     * branch to E = 1.64, which is past Table 5's E ≤ 1 and returns
     * an E' of 1.09, past the [0:1] range Table 5 gives E' as well.
     * Note 5h has signal values exceeding 1.0 during production but
     * does not say what light produces them, and the 109% ceiling
     * itself comes from narrow-range coding (ARIB STD-B67 Table 3,
     * BT.2100-3 Table 9, EBU R 103), not from Table 5. Continuing
     * the log branch above E = 1 is therefore an extrapolation, and
     * an intended one — it is what a camera driving 109% does.
     *
     * Per ITU-R BT.2100-3 Table 5, "Hybrid Log-Gamma (HLG) System"
     */
    hlgOetf(E) {
        E = Math.max(0.0, E);
        if (E <= 1 / 12) {
            // Square-root ("gamma") segment
            return Math.sqrt(3 * E);
        }
        // Logarithmic segment
        return this.a * Math.log(12 * E - this.b) + this.c;
    }

    /**
     * HLG Inverse OETF: Converts non-linear signal value (E') to scene linear light.
     *
     * Per ITU-R BT.2100-3 Table 5, "Hybrid Log-Gamma (HLG) System"
     */
     hlgInverseOetf(E_prime) {
        // Ensure signal is non-negative per BT.2100-3 Table 5
        E_prime = Math.max(0.0, E_prime);
        
        // Apply the piecewise inverse OETF
        if (E_prime <= 0.5) {
            // Square law segment (E ≤ 1/12 in original OETF)
            return (E_prime * E_prime) / 3.0;
        } else {
            // Logarithmic segment (E > 1/12 in original OETF)
            return (Math.exp((E_prime - this.c) / this.a) + this.b) / 12.0;
        }
    }

    /**
     * Calculate the HLG System Gamma according to ITU-R BT.2100-3 and BT.2390.
     * 
     * System gamma is a key parameter that adapts the HLG display response to both:
     * 1. The nominal peak luminance capability (Lw)
     * 2. The viewing environment's surround luminance (Ls)
     * 
     * Two different formulas are specified in ITU-R BT.2100-3 Note 5f.
     * The recommended formula depends on the nominal peak luminance:
     * 
     * - For 400 cd/m² ≤ Lw ≤ 2000 cd/m²: γ = 1.2 + 0.42 log₁₀(Lw/1000)
     * - For Lw < 400 cd/m² or Lw > 2000 cd/m²: γ = 1.2 × κ^(log₂(Lw/1000))
     * 
     * In both cases, additional adjustment for surround luminance is applied:
     * γ = γ × μ^(log₂(Ls/Ls_ref))
     * 
     * Where:
     * - γ_ref = 1.2 (reference gamma at 1000 cd/m² in reference environment)
     * - κ = 1.111 (peak luminance adjustment factor)
     * - μ = 0.98 (surround luminance adjustment factor)
     * - L_ref = 1000 cd/m² (reference nominal peak luminance)
     * - Ls_ref = 5 cd/m² (reference surround luminance)
     * 
     * Note:
     * - Higher nominal peak luminance requires increased gamma (κ > 1)
     * - Brighter viewing environments require decreased gamma (μ < 1)
     */
    calculateSystemGamma(peakLuminance, surroundLuminance) {
        // Select the appropriate formula based on nominal peak luminance
        let baseGamma;
        
        if (peakLuminance >= 400 && peakLuminance <= 2000) {
            // Basic formula for typical production range (400-2000 cd/m²)
            // γ = 1.2 + 0.42 log₁₀(Lw/L_ref)
            baseGamma = this.GAMMA_REF + 0.42 * Math.log10(peakLuminance / this.L_REF);
        } else {
            // Extended formula for displays outside typical range
            // γ = 1.2 × κ^(log₂(Lw/L_ref))
            const kappaFactor = Math.pow(this.KAPPA, Math.log2(peakLuminance / this.L_REF));
            baseGamma = this.GAMMA_REF * kappaFactor;
        }
        
        // Calculate surround luminance adjustment factor (reduces gamma for brighter environments)
        // Following ITU-R BT.2390-12 Section 6.2
        const muFactor = Math.pow(this.MU, Math.log2(surroundLuminance / this.LS_REF));
        
        // Apply surround adjustment to System Gamma
        return baseGamma * muFactor;
    }

    /**
     * Calculate black level lift (β) according to ITU-R BT.2100-3 Table 5.
     * 
     * The black level lift parameter adjusts signal mapping to accommodate non-zero
     * black levels in real displays, ensuring consistent rendering across displays
     * with different black level capabilities.
     * 
     * Per ITU-R BT.2100-3 Table 5, black level lift is defined as:
     *     β = √3(LB/LW)^(1/γ)
     * 
     * Where:
     * - LB is the display luminance for black in cd/m²
     * - LW is nominal peak luminance of the display in cd/m²
     * - γ is the System Gamma
     */
    calculateBlackLevelLift(blackLevel, peakLuminance, gamma) {
        // Early return if black level is 0 or negative (ideal display)
        if (blackLevel <= 0) {
            return 0.0;
        }

        // Calculate black level lift
        return Math.sqrt(3 * Math.pow(blackLevel / peakLuminance, 1 / gamma));
    }

    /**
     * HLG Electro-Optical Transfer Function (EOTF) with black level lift support.
     * 
     * The HLG Reference EOTF consists of these conceptual steps:
     * 1. Apply black level lift if required: E' → max(0,(1-β)E'+β)
     * 2. Inverse OETF: Convert non-linear R', G', B' to scene linear R, G, B
     * 3. Apply OOTF: F_D = α·Y_S^(γ-1)·E, where:
     *    - Y_S is the scene luminance using BT.2020 coefficients
     *    - γ is the variable System Gamma
     *    - α is the user gain in cd/m², representing Lw, so F_D is the
     *      display luminance in cd/m²
     * 
     * From ITU-R BT.2100-3 Table 5. Returns [R_D, G_D, B_D] in cd/m².
     */
    hlgEotf(r_prime, g_prime, b_prime, gamma, peakLuminance, blackLevel = 0, beta = null) {
        // Step 1: Apply black level lift if required (blackLevel > 0)
        const blackLevelLiftEnabled = (blackLevel > 0) || (beta !== null && beta > 0);

        let r_lifted, g_lifted, b_lifted;

        if (blackLevelLiftEnabled) {
            // Calculate beta if not provided
            if (beta === null) {
                beta = this.calculateBlackLevelLift(blackLevel, peakLuminance, gamma);
            }

            // Apply black level lift to each component per ITU-R BT.2100-3 Table 5
            r_lifted = Math.max(0, (1 - beta) * r_prime + beta);
            g_lifted = Math.max(0, (1 - beta) * g_prime + beta);
            b_lifted = Math.max(0, (1 - beta) * b_prime + beta);
        } else {
            // No black level lift
            r_lifted = r_prime;
            g_lifted = g_prime;
            b_lifted = b_prime;
        }

        // Step 2: Inverse OETF - Convert to scene linear
        let r_s = this.hlgInverseOetf(r_lifted);
        let g_s = this.hlgInverseOetf(g_lifted);
        let b_s = this.hlgInverseOetf(b_lifted);

        // Step 3: Calculate scene luminance Y_S using BT.2020 coefficients
        let Y_S = this.Y_R * r_s + this.Y_G * g_s + this.Y_B * b_s;

        // Step 4: Apply the OOTF with α = Lw: F_D = α·Y_S^(γ-1)·E.
        // The single luminance factor preserves colour ratios.
        let ootf_factor = (Y_S > 0) ? Math.pow(Y_S, gamma - 1) : 0;

        // Display luminance per component in cd/m²
        return [peakLuminance * ootf_factor * r_s, peakLuminance * ootf_factor * g_s, peakLuminance * ootf_factor * b_s];
    }

    /**
     * Calculate the HDR Reference White (typically at 75% signal in HLG).
     * 
     * Based on ITU-R BT.2408-9 Section 2.1, HDR Reference White is defined as:
     * "the nominal signal level obtained from an HDR camera and a 100% reflectance 
     * white card resulting in a nominal luminance of 203 cd/m² on a PQ display or 
     * on an HLG display that has a nominal peak luminance capability of 1,000 cd/m²."
     * 
     * For the 1000 cd/m² reference display with γ = 1.2 this function returns
     * 203.15 cd/m², which rounds to the 203 cd/m² of ITU-R BT.2408-9 Table 1. The
     * 75% signal is BT.2408's rounded nominal; 203 cd/m² exactly is 74.99%.
     */
    calculateReferenceWhite(peakLuminance, systemGamma, blackLevel = 0) {
        const signal75 = 0.75; // 75% signal level for HLG Reference White
        const beta = blackLevel > 0 ? this.calculateBlackLevelLift(blackLevel, peakLuminance, systemGamma) : 0;

        // The EOTF returns display luminance in cd/m² (α = Lw)
        const [ref_white] = this.hlgEotf(signal75, signal75, signal75, systemGamma, peakLuminance, blackLevel, beta);

        return ref_white;
    }

    /**
     * Extended Range: Set 100% to peak luminance and calculate the corresponding 109% extended range value.
     */
    computeExtendedRangeValues(peakLuminance, surroundLuminance, blackLevel = 0) {
        const systemGamma = this.calculateSystemGamma(peakLuminance, surroundLuminance);
        const beta = blackLevel > 0 ? this.calculateBlackLevelLift(blackLevel, peakLuminance, systemGamma) : 0;
    
        // Use the precise ITU-R BT.2100-3 109% signal level; the EOTF
        // returns display luminance in cd/m² (α = Lw)
        const [peak109] = this.hlgEotf(this.SIGNAL_LEVEL_109, this.SIGNAL_LEVEL_109, this.SIGNAL_LEVEL_109,
            systemGamma, peakLuminance, blackLevel, beta);
        const refWhite = this.calculateReferenceWhite(peakLuminance, systemGamma, blackLevel);
    
        return {
            lw: peakLuminance,
            gamma: systemGamma,
            peak109: peak109,
            referenceWhite: refWhite,
            blackLevelLift: beta
        };
    }

    /**
     * Solve the nominal peak value Lw from the extended peak value (the
     * luminance the 109% signal is mapped to), so that the full extended
     * range (0-109%) fits within the display peak.
     *
     * This needs an iterative solver because the relationship between the
     * nominal peak luminance and the luminance at 109% cannot be expressed
     * as a direct formula, due to:
     *
     * 1. System gamma depends on the nominal peak luminance (BT.2100-3 Note 5f)
     * 2. The EOTF output at 109% depends on the system gamma
     * 3. The luminance at 109% is the product of nominal peak and EOTF output
     *
     * Lowering the effective nominal peak this way is the "contrast control"
     * adjustment described in BT.2100-3 Note 5f.
     *
     * Note 5f's switch between its two gamma formulas makes γ(Lw)
     * discontinuous at 400 and 2000 cd/m², so a small window of extended
     * peak values (≈3855-3869 at reference surround and zero black level)
     * has no exact solution and another (≈667-671) has two. Bisection relies only on the 109%
     * luminance increasing with Lw within each branch and converges to the
     * largest Lw that fits; a second pass from the Lw = 400 boundary
     * covers the two-solution window.
     */
    solveNominalFromExtendedPeak(extendedPeak, surroundLuminance, blackLevel = 0) {
        // The predicate only needs the 109% luminance; skip computeExtendedRangeValues()'s
        // Reference White computation (a second full EOTF) inside the loop.
        const computePeak109 = lw => {
            const gamma = this.calculateSystemGamma(lw, surroundLuminance);
            const beta = blackLevel > 0 ? this.calculateBlackLevelLift(blackLevel, lw, gamma) : 0;
            const [peak109] = this.hlgEotf(this.SIGNAL_LEVEL_109, this.SIGNAL_LEVEL_109,
                this.SIGNAL_LEVEL_109, gamma, lw, blackLevel, beta);
            return peak109;
        };

        let iterations = 0;
        const largestFittingLw = (lo, hi) => {
            while ((hi - lo) / hi > 1e-7 && iterations < 100) {
                iterations++;
                const mid = (lo + hi) / 2;
                if (computePeak109(mid) <= extendedPeak) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return lo;
        };

        // EOTF[1.0902] ratio is > 1 and < 3 for any reachable gamma, so the
        // solution always lies between extendedPeak / 3 and extendedPeak.
        let lw = largestFittingLw(extendedPeak / 3, extendedPeak);

        // γ steps down across Lw = 400 (Note 5f formula switch), so a
        // fitting Lw just below 400 can coexist with a larger one at or
        // above it; prefer the larger so mode swaps round-trip.
        if (lw < 400 && computePeak109(400) <= extendedPeak) {
            lw = largestFittingLw(400, extendedPeak);
        }

        return this.computeExtendedRangeValues(lw, surroundLuminance, blackLevel);
    }
}

//=============================================================================
// UI CONTROLLER
//=============================================================================

/**
 * Manages UI elements and user interaction
 */
class UIController {
    constructor(calculator) {
        this.calculator = calculator;

        // Set by the application after construction; in extended-peak
        // mode that controller owns the derived table cells
        this.extendedRangeController = null;

        // Slider elements
        this.peakLuminanceSlider = document.getElementById("peakLuminance");
        this.surroundLuminanceSlider = document.getElementById("surroundLuminance");
        this.blackLevelSlider = document.getElementById("blackLevel");

        // Input field elements
        this.peakLuminanceInput = document.getElementById("peakLuminanceInput");
        this.surroundLuminanceInput = document.getElementById("surroundLuminanceInput");
        this.blackLevelInput = document.getElementById("blackLevelInput");

        // Display value elements
        this.peakLuminanceValue = document.getElementById("peakLuminanceValue");
        this.surroundLuminanceValue = document.getElementById("surroundLuminanceValue");
        this.blackLevelValue = document.getElementById("blackLevelValue");
        this.blackLiftValue = document.getElementById("blackLiftValue");

        this.systemGammaValue = document.getElementById("systemGamma");
        this.referenceWhiteValue = document.getElementById("referenceWhite");

        // Track last values to avoid unnecessary updates
        this.lastPeakLuminance = null;
        this.lastSurroundLuminance = null;
        this.lastBlackLevel = null;
        this.blackLevelEnabled = false;
    }

    /**
     * Set up event listeners for UI elements
     */
    setupEventListeners() {
        // Peak Luminance slider and input events
        this.peakLuminanceSlider.addEventListener("input", () => {
            this.syncInputFromSlider('peak');
            this.updateDisplay(true);
            this.emitUpdateEvent();
        });

        this.peakLuminanceInput.addEventListener("input", () => {
            this.syncSliderFromInput('peak');
            this.updateDisplay(true);
            this.emitUpdateEvent();
        });

        this.peakLuminanceInput.addEventListener("blur", () => {
            this.validateInput('peak');
        });

        // Surround Luminance slider and input events
        this.surroundLuminanceSlider.addEventListener("input", () => {
            this.syncInputFromSlider('surround');
            this.updateDisplay(true);
            this.emitUpdateEvent();
        });

        this.surroundLuminanceInput.addEventListener("input", () => {
            this.syncSliderFromInput('surround');
            this.updateDisplay(true);
            this.emitUpdateEvent();
        });

        this.surroundLuminanceInput.addEventListener("blur", () => {
            this.validateInput('surround');
        });

        // Black Level slider and input events
        this.blackLevelSlider.addEventListener("input", () => {
            this.syncInputFromSlider('blackLevel');
            this.updateDisplay(true);
            this.emitUpdateEvent();
        });

        this.blackLevelInput.addEventListener("input", () => {
            this.syncSliderFromInput('blackLevel');
            this.updateDisplay(true);
            this.emitUpdateEvent();
        });

        this.blackLevelInput.addEventListener("blur", () => {
            this.validateInput('blackLevel');
        });
    }

    /**
     * Sync input field from slider value
     */
    syncInputFromSlider(type) {
        if (type === 'peak') {
            this.peakLuminanceInput.value = this.peakLuminanceSlider.value;
        } else if (type === 'surround') {
            this.surroundLuminanceInput.value = this.surroundLuminanceSlider.value;
        } else if (type === 'blackLevel') {
            this.blackLevelInput.value = this.blackLevelSlider.value;
        }
    }

    /**
     * Sync slider from input field value
     */
    syncSliderFromInput(type) {
        if (type === 'peak') {
            const value = parseFloat(this.peakLuminanceInput.value);
            if (!isNaN(value)) {
                this.peakLuminanceSlider.value = Math.max(
                    this.peakLuminanceSlider.min, 
                    Math.min(this.peakLuminanceSlider.max, value)
                );
            }
        } else if (type === 'surround') {
            const value = parseFloat(this.surroundLuminanceInput.value);
            if (!isNaN(value)) {
                this.surroundLuminanceSlider.value = Math.max(
                    this.surroundLuminanceSlider.min, 
                    Math.min(this.surroundLuminanceSlider.max, value)
                );
            }
        } else if (type === 'blackLevel') {
            const value = parseFloat(this.blackLevelInput.value);
            if (!isNaN(value)) {
                this.blackLevelSlider.value = Math.max(
                    this.blackLevelSlider.min, 
                    Math.min(this.blackLevelSlider.max, value)
                );
            }
        }
    }

    /**
     * Validate and constrain input values
     */
    validateInput(type) {
        if (type === 'peak') {
            const min = parseFloat(this.peakLuminanceSlider.min);
            const max = parseFloat(this.peakLuminanceSlider.max);
            const value = parseFloat(this.peakLuminanceInput.value);
            
            if (isNaN(value) || value < min || value > max) {
                this.peakLuminanceInput.value = this.peakLuminanceSlider.value;
            }
        } else if (type === 'surround') {
            const min = parseFloat(this.surroundLuminanceSlider.min);
            const max = parseFloat(this.surroundLuminanceSlider.max);
            const value = parseFloat(this.surroundLuminanceInput.value);
            
            if (isNaN(value) || value < min || value > max) {
                this.surroundLuminanceInput.value = this.surroundLuminanceSlider.value;
            }
        } else if (type === 'blackLevel') {
            const min = parseFloat(this.blackLevelSlider.min);
            const max = parseFloat(this.blackLevelSlider.max);
            const value = parseFloat(this.blackLevelInput.value);
            
            if (isNaN(value) || value < min || value > max) {
                this.blackLevelInput.value = this.blackLevelSlider.value;
            }
        }
    }

    /**
     * Enable or disable black level
     */
    setBlackLevelEnabled(enabled) {
        this.blackLevelEnabled = enabled;

        // Show/hide black level UI elements
        const blackLevelSliderContainer = document.getElementById("blackLevelSliderContainer");
        const blackLevelRow = document.getElementById("blackLevelRow");
        const blackLiftRow = document.getElementById("blackLiftRow");
        const blackLevelFormulaSection = document.getElementById("blackLevelFormulaSection");

        if (enabled) {
            blackLevelSliderContainer.style.display = "block";
            blackLevelRow.style.display = "table-row";
            blackLiftRow.style.display = "table-row";
            blackLevelFormulaSection.style.display = "block";
        } else {
            blackLevelSliderContainer.style.display = "none";
            blackLevelRow.style.display = "none";
            blackLiftRow.style.display = "none";
            blackLevelFormulaSection.style.display = "none";
        }

        // Update display with new state
        this.updateDisplay(true);
        this.emitUpdateEvent();
    }

    /**
     * Update display with current values
     * 
     * @param {boolean} forceUpdate - Whether to force update even if values haven't changed
     * @returns {Object|null} Calculated values, or null when unchanged or deferred
     */
    updateDisplay(forceUpdate = false) {
        const peakLuminance = parseFloat(this.peakLuminanceSlider.value);
        const surroundLuminance = parseFloat(this.surroundLuminanceSlider.value);
        const blackLevel = this.blackLevelEnabled ? parseFloat(this.blackLevelSlider.value) : 0;

        // Skip update if values haven't changed (unless forced)
        if (!forceUpdate &&
            peakLuminance === this.lastPeakLuminance &&
            surroundLuminance === this.lastSurroundLuminance &&
            blackLevel === this.lastBlackLevel) {
            return null;
        }

        this.lastPeakLuminance = peakLuminance;
        this.lastSurroundLuminance = surroundLuminance;
        this.lastBlackLevel = blackLevel;

        // Update the input cells
        this.surroundLuminanceValue.textContent = `${surroundLuminance} cd/m²`;
        if (this.blackLevelEnabled) {
            this.blackLevelValue.textContent = `${blackLevel.toFixed(4)} cd/m²`;
        }

        // In extended-peak mode the ExtendedRangeController owns the
        // peak, System Gamma, Reference White and black lift cells:
        // the slider value is the extended peak, not the nominal Lw
        // these would be derived from.
        if (this.extendedRangeController && this.extendedRangeController.sliderSetsExtendedPeak) {
            return null;
        }

        // Calculate the standard values
        const systemGamma = this.calculator.calculateSystemGamma(peakLuminance, surroundLuminance);
        const blackLevelLift = blackLevel > 0 ? this.calculator.calculateBlackLevelLift(blackLevel, peakLuminance, systemGamma) : 0;
        const referenceWhite = this.calculator.calculateReferenceWhite(peakLuminance, systemGamma, blackLevel);

        // Update the display
        this.peakLuminanceValue.textContent = `${peakLuminance} cd/m²`;
        this.systemGammaValue.textContent = systemGamma.toFixed(2);
        this.referenceWhiteValue.textContent = `${referenceWhite.toFixed(0)} cd/m²`;
        if (this.blackLevelEnabled) {
            this.blackLiftValue.textContent = blackLevelLift.toFixed(6);
        }

        // Return calculated values for other components to use
        return {
            peakLuminance,
            surroundLuminance,
            blackLevel,
            systemGamma,
            blackLevelLift,
            referenceWhite
        };
    }

    /**
     * Emit an event when values are updated
     */
    emitUpdateEvent() {
        const peakLuminance = parseFloat(this.peakLuminanceSlider.value);
        const surroundLuminance = parseFloat(this.surroundLuminanceSlider.value);
        const blackLevel = this.blackLevelEnabled ? parseFloat(this.blackLevelSlider.value) : 0;

        // Create and dispatch the event
        const event = new CustomEvent('hlg-values-updated', {
            detail: {
                peakLuminance,
                surroundLuminance,
                blackLevel,
                blackLevelEnabled: this.blackLevelEnabled
            }
        });

        document.dispatchEvent(event);
    }

    /**
     * Get current values
     * 
     * @returns {Object} Current parameter values
     */
    getCurrentValues() {
        return {
            peakLuminance: parseFloat(this.peakLuminanceSlider.value),
            surroundLuminance: parseFloat(this.surroundLuminanceSlider.value),
            blackLevel: this.blackLevelEnabled ? parseFloat(this.blackLevelSlider.value) : 0,
            blackLevelEnabled: this.blackLevelEnabled
        };
    }

    /**
     * Update display with provided values (for extended range mode)
     * 
     * @param {number} gamma - System gamma value to display
     * @param {number} referenceWhite - Reference white value to display
     * @param {number} blackLevelLift - Black level lift value to display
     */
    displayDerivedValues(gamma, referenceWhite, blackLevelLift = 0) {
        this.systemGammaValue.textContent = gamma.toFixed(2);
        this.referenceWhiteValue.textContent = `${referenceWhite.toFixed(0)} cd/m²`;

        if (this.blackLevelEnabled) {
            this.blackLiftValue.textContent = blackLevelLift.toFixed(6);
        }
    }
}

//=============================================================================
// BLACK LEVEL CONTROLLER
//=============================================================================

/**
 * Manages black level functionality
 */
class BlackLevelController {
    constructor(uiController) {
        this.uiController = uiController;

        // Buttons
        this.blackLevelOffBtn = document.getElementById("blackLevelOffBtn");
        this.blackLevelOnBtn = document.getElementById("blackLevelOnBtn");

        // Black level slider
        this.blackLevelSlider = document.getElementById("blackLevel");

        // Black level buttons container
        this.blackLevelButtons = document.getElementById("blackLevelButtons");

        // Preset values
        this.BLACK_LEVEL_PRESETS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
    }

    /**
     * Initialize black level controller
     */
    initialize() {
        this.setupButtons();
        this.setupBlackLevelPresets();
    }

    /**
     * Set up black level toggle buttons
     */
    setupButtons() {
        this.blackLevelOffBtn.addEventListener("click", () => {
            this.setBlackLevelEnabled(false);
        });

        this.blackLevelOnBtn.addEventListener("click", () => {
            this.setBlackLevelEnabled(true);
        });
    }

    /**
     * Set up black level preset buttons
     */
    setupBlackLevelPresets() {
        this.blackLevelButtons.innerHTML = "";

        this.BLACK_LEVEL_PRESETS.forEach(value => {
            const button = document.createElement("button");
            button.textContent = value.toFixed(value < 0.01 ? 3 : 2);

            button.onclick = () => {
                this.blackLevelSlider.value = value;
                this.blackLevelSlider.dispatchEvent(new Event("input"));
                this.updateButtonHighlight(value);
            };

            this.blackLevelButtons.appendChild(button);
        });

        // Set initial highlight
        this.updateButtonHighlight(parseFloat(this.blackLevelSlider.value));

        // Update highlight when slider changes
        this.blackLevelSlider.addEventListener("input", () => {
            this.updateButtonHighlight(parseFloat(this.blackLevelSlider.value));
        });
    }

    /**
     * Set black level enabled/disabled
     */
    setBlackLevelEnabled(enabled) {
        // Toggle button active states
        if (enabled) {
            this.blackLevelOffBtn.classList.remove("active");
            this.blackLevelOnBtn.classList.add("active");
        } else {
            this.blackLevelOffBtn.classList.add("active");
            this.blackLevelOnBtn.classList.remove("active");
        }

        // Update UI controller
        this.uiController.setBlackLevelEnabled(enabled);
    }

    /**
     * Update button highlighting based on slider value
     */
    updateButtonHighlight(value) {
        const buttons = this.blackLevelButtons.querySelectorAll("button");

        buttons.forEach(button => {
            if (Math.abs(parseFloat(button.textContent) - value) < 0.0001) {
                button.classList.add("active");
            } else {
                button.classList.remove("active");
            }
        });
    }
}

//=============================================================================
// FORMULA RENDERER
//=============================================================================

/**
 * Handles rendering of mathematical formulas using KaTeX
 */
class FormulaRenderer {
    constructor(calculator) {
        this.calculator = calculator;

        // Formula display elements
        // Set by the application after construction; when Extended Range
        // is active that controller owns the formula rendering
        this.extendedRangeController = null;

        this.mathGamma = document.getElementById("mathGamma");
        this.mathRefWhite = document.getElementById("mathRefWhite");
        this.mathBlackLift = document.getElementById("mathBlackLift");
        this.mathExtendedRange = document.getElementById("mathExtendedRange");
    }

    /**
     * Initialize KaTeX rendering
     */
    initialize() {
        this.setupEventListeners();
        this.renderStaticFormulas();
    }

    /**
     * Set up event listeners for formula updates
     */
    setupEventListeners() {
        document.addEventListener('hlg-values-updated', (event) => {
            // In Extended Range mode the ExtendedRangeController re-renders
            // the formulas with the extended data (and, in extended-peak
            // mode, the solved nominal peak), so this render would be
            // overwritten unread — and in extended-peak mode the raw
            // slider value is not the nominal Lw this render assumes.
            if (this.extendedRangeController && this.extendedRangeController.extendedRangeActive) {
                return;
            }
            this.updateFormulas(
                event.detail.peakLuminance,
                event.detail.surroundLuminance,
                event.detail.blackLevel,
                event.detail.blackLevelEnabled
            );
        });
    }

    /**
     * Render all static formulas in the document
     */
    renderStaticFormulas() {
        // Find and render all static formulas
        document.querySelectorAll(".formula-container:not(#mathGamma):not(#mathRefWhite):not(#mathBlackLift)").forEach(el => {
            if (window.katex) {
                katex.render(el.textContent, el, { throwOnError: false });
            }
        });
    }

    /**
     * Update dynamic formulas based on current values
     * 
     * @param {number} peakLuminance - Current peak luminance
     * @param {number} surroundLuminance - Current surround luminance
     * @param {number} blackLevel - Current black level
     * @param {boolean} blackLevelEnabled - Whether black level is enabled
     * @param {Object} extendedRangeData - Extended range calculation results (optional)
     */
    updateFormulas(peakLuminance, surroundLuminance, blackLevel = 0, blackLevelEnabled = false, extendedRangeData = null) {
        if (!window.katex) return;

        let peak100 = peakLuminance;
        let systemGamma, referenceWhite, blackLevelLift = 0;

        // Use extended range data if provided, otherwise calculate standard values
        if (extendedRangeData) {
            peak100 = extendedRangeData.lw;
            systemGamma = extendedRangeData.gamma;
            referenceWhite = extendedRangeData.referenceWhite;
            blackLevelLift = extendedRangeData.blackLevelLift || 0;
        } else {
            systemGamma = this.calculator.calculateSystemGamma(peak100, surroundLuminance);

            if (blackLevelEnabled) {
                blackLevelLift = this.calculator.calculateBlackLevelLift(blackLevel, peak100, systemGamma);
            }

            referenceWhite = this.calculator.calculateReferenceWhite(peak100, systemGamma, blackLevelEnabled ? blackLevel : 0);
        }

        // Determine which gamma formula is used based on peak luminance
        let gammaFormula;
        let usesExtendedGammaFormula = (peak100 < 400 || peak100 > 2000);
        
        if (usesExtendedGammaFormula) {
            // Extended range formula (for Lw < 400 or Lw > 2000)
            gammaFormula = String.raw`
                \gamma = 1.2 \times 1.111^{\log_2\left(\frac{${peak100.toFixed(2)}}{1000}\right)} 
                \times 0.98^{\log_2\left(\frac{${surroundLuminance.toFixed(2)}}{5}\right)} = ${systemGamma.toFixed(3)}
            `;
        } else {
            // Basic formula (for 400 ≤ Lw ≤ 2000)
            gammaFormula = String.raw`
                \gamma = \left(1.2 + 0.42 \log_{10}\left(\frac{${peak100.toFixed(2)}}{1000}\right)\right) 
                \times 0.98^{\log_2\left(\frac{${surroundLuminance.toFixed(2)}}{5}\right)} = ${systemGamma.toFixed(3)}
            `;
        }

        // Update Reference White formula
        let refWhiteFormula = String.raw`
            L_{\text{ref white}} = ${peak100.toFixed(2)} \times \left( \frac{\exp\left(\frac{0.75-0.5599}{0.1788}\right)+0.2847}{12} \right)^{${systemGamma.toFixed(3)}} = ${referenceWhite.toFixed(2)} \text{ cd/m}^{2}
        `;

        // Add black level lift formula if enabled
        if (blackLevelEnabled) {
            const blackLiftFormula = String.raw`
                \beta = \sqrt{3 \left(\frac{${blackLevel.toFixed(4)}}{${peak100.toFixed(2)}}\right)^{1/${systemGamma.toFixed(3)}}} = ${blackLevelLift.toFixed(6)}
            `;
        
            katex.render(blackLiftFormula, this.mathBlackLift, { throwOnError: false });
        
            // Add black level lift to Reference White formula. The
            // lift appears as the symbol beta (its value is derived
            // in the formula directly above), as in the BT.2100-3
            // Table 5 EOTF; substituting the six-decimal value twice
            // would make the term wider than a phone screen.
            refWhiteFormula = String.raw`
                L_{\text{ref white}} = ${peak100.toFixed(2)} \times \left( \text{OETF}^{-1}\left[ \max(0, (1-\beta) \times 0.75 + \beta) \right] \right)^{${systemGamma.toFixed(3)}} = ${referenceWhite.toFixed(2)} \text{ cd/m}^{2}
            `;
        }

        // Add Extended Range formula when extendedRangeData is provided
        if (extendedRangeData) {
            let extendedRangeFormula;

            if (blackLevelEnabled && extendedRangeData.blackLevelLift > 0) {
                // beta symbolically, as in the Reference White formula
                extendedRangeFormula = String.raw`
                    L_{109\%} = ${peak100.toFixed(2)} \times \left( \text{OETF}^{-1}\left[ \max(0, (1-\beta) \times ${this.calculator.SIGNAL_LEVEL_109.toFixed(6)} + \beta) \right] \right)^{${systemGamma.toFixed(3)}} = ${extendedRangeData.peak109.toFixed(2)} \text{ cd/m}^{2}
                `;
            } else {
                extendedRangeFormula = String.raw`
                    L_{109\%} = ${peak100.toFixed(2)} \times \left( \text{OETF}^{-1}\left[ ${this.calculator.SIGNAL_LEVEL_109.toFixed(6)} \right] \right)^{${systemGamma.toFixed(3)}} = ${peak100.toFixed(2)} \times ${(extendedRangeData.peak109/peak100).toFixed(6)} = ${extendedRangeData.peak109.toFixed(2)} \text{ cd/m}^{2}
                `;
            }

            katex.render(extendedRangeFormula, this.mathExtendedRange, { throwOnError: false });
        }

        // Render using KaTeX
        katex.render(gammaFormula, this.mathGamma, { throwOnError: false });
        katex.render(refWhiteFormula, this.mathRefWhite, { throwOnError: false });
    }
}

//=============================================================================
// PRESET CONTROLLER
//=============================================================================

/**
 * Manages preset buttons for peak and surround luminance
 */
class PresetController {
    constructor(uiController) {
        this.uiController = uiController;

        // Preset values
        this.PEAK_LUMINANCE_PRESETS = [400, 500, 600, 1000, 2000, 3000, 4000, 5000, 10000];
        this.SURROUND_LUMINANCE_PRESETS = [0.05, 1, 5, 10, 20, 50, 100, 200, 500];

        // Button containers
        this.peakButtonsContainer = document.getElementById("peakLuminanceButtons");
        this.surroundButtonsContainer = document.getElementById("surroundLuminanceButtons");

        // Slider elements
        this.peakLuminanceSlider = document.getElementById("peakLuminance");
        this.surroundLuminanceSlider = document.getElementById("surroundLuminance");
    }

    /**
     * Initialize preset buttons
     */
    initialize() {
        this.createPresetButtons();
        this.setupEventListeners();
    }

    /**
     * Create preset buttons for both sliders
     */
    createPresetButtons() {
        this.createButtonsForContainer(
            this.peakButtonsContainer,
            this.PEAK_LUMINANCE_PRESETS,
            this.peakLuminanceSlider
        );

        this.createButtonsForContainer(
            this.surroundButtonsContainer,
            this.SURROUND_LUMINANCE_PRESETS,
            this.surroundLuminanceSlider
        );

        // Highlight the initial active buttons
        this.updateButtonHighlight(
            this.peakButtonsContainer,
            this.PEAK_LUMINANCE_PRESETS,
            parseFloat(this.peakLuminanceSlider.value)
        );

        this.updateButtonHighlight(
            this.surroundButtonsContainer,
            this.SURROUND_LUMINANCE_PRESETS,
            parseFloat(this.surroundLuminanceSlider.value)
        );
    }

    /**
     * Create preset buttons for a specific container
     * 
     * @param {HTMLElement} container - Button container element
     * @param {Array<number>} values - Preset values
     * @param {HTMLInputElement} slider - Associated slider element
     */
    createButtonsForContainer(container, values, slider) {
        container.innerHTML = "";

        values.forEach(value => {
            const button = document.createElement("button");
            button.textContent = value;

            button.onclick = () => {
                slider.value = value;
                slider.dispatchEvent(new Event("input"));
                this.updateButtonHighlight(container, values, value);
            };

            container.appendChild(button);
        });
    }

    /**
     * Set up event listeners for slider changes
     */
    setupEventListeners() {
        this.peakLuminanceSlider.addEventListener("input", () => {
            this.updateButtonHighlight(
                this.peakButtonsContainer,
                this.PEAK_LUMINANCE_PRESETS,
                parseFloat(this.peakLuminanceSlider.value)
            );
        });

        this.surroundLuminanceSlider.addEventListener("input", () => {
            this.updateButtonHighlight(
                this.surroundButtonsContainer,
                this.SURROUND_LUMINANCE_PRESETS,
                parseFloat(this.surroundLuminanceSlider.value)
            );
        });
    }

    /**
     * Update button highlighting based on current slider value
     * 
     * @param {HTMLElement} container - Button container
     * @param {Array<number>} values - Preset values
     * @param {number} currentValue - Current slider value
     */
    updateButtonHighlight(container, values, currentValue) {
        const buttons = container.querySelectorAll("button");

        buttons.forEach(button => {
            if (parseFloat(button.textContent) === parseFloat(currentValue)) {
                button.classList.add("active");
            } else {
                button.classList.remove("active");
            }
        });
    }
}

//=============================================================================
// SLIDER RESET CONTROLLER
//=============================================================================

/**
 * The reset beside each display slider's value box. Each puts back the
 * one value it stands beside and leaves the rest of the display where
 * the reader put it: the reference display's nominal peak L_W and its
 * surround (BT.2390-12 section 6.2), and the black level this page
 * opens at, which no standard fixes.
 *
 * The graph views repeat these sliders and their resets click these,
 * so a value goes back everywhere it is shown at once.
 */
class SliderResetController {
    constructor(calculator, uiController, extendedRangeController) {
        this.calculator = calculator;
        this.uiController = uiController;
        this.extendedRangeController = extendedRangeController;

        this.BLACK_LEVEL_DEFAULT = 0.005;

        this.peakSlider = document.getElementById("peakLuminance");
        this.surroundSlider = document.getElementById("surroundLuminance");
        this.blackLevelSlider = document.getElementById("blackLevel");

        this.peakButton = document.getElementById("peakLuminanceResetBtn");
        this.surroundButton = document.getElementById("surroundLuminanceResetBtn");
        this.blackLevelButton = document.getElementById("blackLevelResetBtn");
    }

    initialize() {
        this.peakButton.addEventListener("click",
            () => this.setSlider(this.peakSlider, this.referencePeak()));
        this.surroundButton.addEventListener("click",
            () => this.setSlider(this.surroundSlider, this.calculator.LS_REF));
        this.blackLevelButton.addEventListener("click",
            () => this.setSlider(this.blackLevelSlider, this.BLACK_LEVEL_DEFAULT));
    }

    /**
     * The reference display's peak in whichever quantity the slider is
     * currently setting: its nominal L_W, or in extended-peak mode the
     * luminance that display maps the 109% signal to, which the
     * surround and black level place (ExtendedRangeController converts
     * between the two the same way).
     */
    referencePeak() {
        if (!this.extendedRangeController.sliderSetsExtendedPeak) return this.calculator.L_REF;
        const { surroundLuminance, blackLevel } = this.uiController.getCurrentValues();
        return Math.round(this.calculator.computeExtendedRangeValues(
            this.calculator.L_REF, surroundLuminance, blackLevel).peak109);
    }

    /**
     * Set a value on the slider that owns it, clamped to its range as
     * the calculator's own value boxes clamp
     */
    setSlider(slider, value) {
        const min = parseFloat(slider.min);
        const max = parseFloat(slider.max);
        slider.value = Math.min(max, Math.max(min, value));
        slider.dispatchEvent(new Event("input"));
    }
}

//=============================================================================
// EXTENDED RANGE CONTROLLER
//=============================================================================

/**
 * Manages extended range (0-109%) functionality
 */
class ExtendedRangeController {
    constructor(calculator, uiController, formulaRenderer) {
        this.calculator = calculator;
        this.uiController = uiController;
        this.formulaRenderer = formulaRenderer;
    
        // State
        this.extendedRangeActive = false;
        this.sliderSetsExtendedPeak = false;

        // Elements
        this.nominalRangeButton = document.getElementById("nominalRangeButton");
        this.extendedRangeButton = document.getElementById("extendedRangeButton");
        this.extendedPeakRow = document.getElementById("extendedPeakRow");
        this.extendedPeakValue = document.getElementById("extendedPeakValue");
        this.extendedRangeFormulaSection = document.getElementById("extendedRangeFormulaSection");
        this.extendedRangeNoteSection = document.getElementById("extendedRangeNoteSection");

        // Nominal/extended peak toggle elements
        this.peakNominalOrExtendedToggle = document.getElementById("peakNominalOrExtendedToggle");
        this.peakNominalBtn = document.getElementById("peakNominalBtn");
        this.peakExtendedBtn = document.getElementById("peakExtendedBtn");
        this.peakLuminanceTitle = document.getElementById("peakLuminanceTitle");
        this.peakLuminanceSignal = document.getElementById("peakLuminanceSignal");
        this.nominalPeakRow = document.getElementById("nominalPeakRow");
        this.nominalPeakRoleTag = document.getElementById("nominalPeakRole");
        this.extendedPeakRoleTag = document.getElementById("extendedPeakRole");
        this.peakLuminanceValueCell = document.getElementById("peakLuminanceValue");
        this.extendedPeakNoteSection = document.getElementById("extendedPeakNoteSection");
        this.extendedRangeDetails = document.getElementById("extendedRangeDetails");
        this.extendedRangeSummaryText = document.getElementById("extendedRangeSummaryText");
        this.peakSlider = document.getElementById("peakLuminance");
        this.peakInput = document.getElementById("peakLuminanceInput");
    }

    /**
     * Initialize extended range controller
     */
    initialize() {
        this.setupEventListeners();
        this.hideExtendedRangeElements();
    }

    /**
     * Set up event listeners for extended range buttons
     */
    setupEventListeners() {
        // Toggle buttons
        this.nominalRangeButton.addEventListener("click", () => this.setNominalRangeMode());
        this.extendedRangeButton.addEventListener("click", () => this.setExtendedRangeMode());

        // Nominal/extended peak toggle buttons
        this.peakNominalBtn.addEventListener("click", () => this.setPeakMode(false));
        this.peakExtendedBtn.addEventListener("click", () => this.setPeakMode(true));

        // Update values when parameters change
        document.addEventListener('hlg-values-updated', (event) => {
            if (this.extendedRangeActive) {
                if (this.sliderSetsExtendedPeak) {
                    this.updateExtendedPeakValues(
                        event.detail.peakLuminance,
                        event.detail.surroundLuminance,
                        event.detail.blackLevel,
                        event.detail.blackLevelEnabled
                    );
                } else {
                    this.updateExtendedRangeValues(
                        event.detail.peakLuminance,
                        event.detail.surroundLuminance,
                        event.detail.blackLevel,
                        event.detail.blackLevelEnabled
                    );
                }
            }
        });
    }

    /**
     * Hide extended range UI elements
     */
    hideExtendedRangeElements() {
        this.extendedPeakRow.style.display = "none";
        this.extendedRangeFormulaSection.style.display = "none";
        this.extendedRangeNoteSection.style.display = "none";
        this.extendedPeakNoteSection.style.display = "none";
        this.extendedRangeDetails.style.display = "none";
        this.peakNominalOrExtendedToggle.style.display = "none";
    }

    /**
     * Switch to nominal range mode (0-100%)
     */
    setNominalRangeMode() {
        this.extendedRangeActive = false;
        this.nominalRangeButton.classList.add("active");
        this.extendedRangeButton.classList.remove("active");

        // Leaving Extended Range: the extended-peak mode no longer applies.
        // Revert the slider to nominal, carrying the solved configuration over.
        if (this.sliderSetsExtendedPeak) {
            this.setPeakMode(false);
        }

        this.hideExtendedRangeElements();

        // Update with current values
        const { peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled } = this.uiController.getCurrentValues();
        this.formulaRenderer.updateFormulas(peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled);
    }

    /**
     * Switch to extended range mode (0-109%)
     */
    setExtendedRangeMode() {
        this.extendedRangeActive = true;
        this.extendedRangeButton.classList.add("active");
        this.nominalRangeButton.classList.remove("active");

        // Show extended range elements
        this.extendedPeakRow.style.display = "table-row";
        this.extendedRangeFormulaSection.style.display = "block";
        this.peakNominalOrExtendedToggle.style.display = "flex";
        this.extendedRangeDetails.style.display = "block";
        this.applyPeakModePresentation();

        // Update values with current settings
        const { peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled } = this.uiController.getCurrentValues();
        if (this.sliderSetsExtendedPeak) {
            this.updateExtendedPeakValues(peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled);
        } else {
            this.updateExtendedRangeValues(peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled);
        }
    }

    /**
     * Switch what the peak luminance slider sets: the nominal peak Lw
     * (setsExtendedPeak = false, the default) or the extended peak,
     * i.e. the luminance the 109% signal is mapped to (setsExtendedPeak = true).
     *
     * The swap carries the same display configuration across (like a
     * currency converter): the slider value changes to the equivalent
     * quantity, so the displayed γ / Reference White do not jump.
     *
     * Near the ends of the slider ranges the equivalent value can fall
     * outside the target range (the extended maximum of 25000 reaches
     * past the 109% luminance of the nominal maximum, ≈23200 at
     * reference surround, and the extended minimum of 100 solves to a
     * nominal peak below 100): it is then clamped to the bounds and
     * the configuration follows the clamped value.
     */
    setPeakMode(setsExtendedPeak) {
        if (setsExtendedPeak === this.sliderSetsExtendedPeak) return;

        const { peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled } = this.uiController.getCurrentValues();
        const effectiveBlackLevel = blackLevelEnabled ? blackLevel : 0;

        let newSliderValue;
        if (setsExtendedPeak) {
            // Slider held nominal Lw -> becomes the corresponding 109% luminance
            newSliderValue = this.calculator.computeExtendedRangeValues(peakLuminance, surroundLuminance, effectiveBlackLevel).peak109;
        } else {
            // Slider held the extended peak -> becomes the solved nominal Lw
            newSliderValue = this.calculator.solveNominalFromExtendedPeak(peakLuminance, surroundLuminance, effectiveBlackLevel).lw;
        }

        this.sliderSetsExtendedPeak = setsExtendedPeak;
        this.applyPeakModePresentation();

        // Clamp explicitly (a range input would clamp silently): the
        // slider ranges are fixed, so an equivalent value outside the
        // target range snaps to the nearest bound
        const min = parseFloat(this.peakSlider.min);
        const max = parseFloat(this.peakSlider.max);
        this.peakSlider.value = Math.min(max, Math.max(min, Math.round(newSliderValue)));
        this.peakSlider.dispatchEvent(new Event("input"));
    }

    /**
     * Apply the presentation for the current peak mode: button states,
     * slider label and range, green/blue (set/calculated)
     * colour roles, and which extended-range note is shown.
     */
    applyPeakModePresentation() {
        const setsExtendedPeak = this.sliderSetsExtendedPeak;

        this.peakNominalBtn.classList.toggle("active", !setsExtendedPeak);
        this.peakExtendedBtn.classList.toggle("active", setsExtendedPeak);

        // Slider label and range (109% of the 10000 cd/m² nominal maximum
        // needs headroom above 10000). The quantity heads the label and
        // the signal level it belongs to sits under it, so neither has
        // to share a wrapped line with the unit
        this.peakLuminanceTitle.textContent = setsExtendedPeak
            ? "Extended Peak Luminance (cd/m²)"
            : "Nominal Peak Luminance (cd/m²)";
        this.peakLuminanceSignal.textContent = setsExtendedPeak
            ? "109% signal"
            : "100% signal";
        const max = setsExtendedPeak ? 25000 : 10000;
        this.peakSlider.max = max;
        this.peakInput.max = max;

        // Table: green = set by the slider, blue = calculated. In the
        // extended-peak mode the 109% row is the input and the nominal
        // row is solved, so the colours and the role pills swap;
        // the rows keep their fixed order and labels.
        this.peakLuminanceValueCell.classList.toggle("output-value-slider", !setsExtendedPeak);
        this.peakLuminanceValueCell.classList.toggle("output-value-calculated", setsExtendedPeak);
        this.extendedPeakValue.classList.toggle("output-value-slider", setsExtendedPeak);
        this.extendedPeakValue.classList.toggle("output-value-calculated", !setsExtendedPeak);

        this.nominalPeakRoleTag.textContent = setsExtendedPeak ? "Calculated" : "Input";
        this.nominalPeakRoleTag.classList.toggle("role-input", !setsExtendedPeak);
        this.nominalPeakRoleTag.classList.toggle("role-calculated", setsExtendedPeak);
        this.extendedPeakRoleTag.textContent = setsExtendedPeak ? "Input" : "Calculated";
        this.extendedPeakRoleTag.classList.toggle("role-input", setsExtendedPeak);
        this.extendedPeakRoleTag.classList.toggle("role-calculated", !setsExtendedPeak);

        // Show note about Extended Range
        this.extendedRangeNoteSection.style.display = (this.extendedRangeActive && !setsExtendedPeak) ? "block" : "none";
        this.extendedPeakNoteSection.style.display = (this.extendedRangeActive && setsExtendedPeak) ? "block" : "none";

        // Summary line carries the one-line lesson for the current peak mode;
        // the open/closed state is left alone so an opened note stays open
        // while the peak mode is toggled.
        this.extendedRangeSummaryText.innerHTML = setsExtendedPeak
            ? "The nominal peak L<sub>W</sub> is calculated from the extended peak value"
            : "Extended Range does not change system gamma or HDR Reference White";

    }

    /**
     * Update values when the slider sets the extended peak:
     * solve the nominal Lw that fits the full extended range, then
     * present the solved configuration.
     *
     * @param {number} extendedPeak - Extended peak luminance (slider value)
     * @param {number} surroundLuminance - Current surround luminance
     * @param {number} blackLevel - Current black level
     * @param {boolean} blackLevelEnabled - Whether black level is enabled
     */
    updateExtendedPeakValues(extendedPeak, surroundLuminance, blackLevel = 0, blackLevelEnabled = false) {
        const effectiveBlackLevel = blackLevelEnabled ? blackLevel : 0;

        const data = this.calculator.solveNominalFromExtendedPeak(extendedPeak, surroundLuminance, effectiveBlackLevel);

        // Update the main table values (the nominal row is now solved)
        this.uiController.displayDerivedValues(
            data.gamma,
            data.referenceWhite,
            data.blackLevelLift
        );
        this.peakLuminanceValueCell.textContent = `${data.lw.toFixed(0)} cd/m²`;
        // Show the solved configuration's actual 109% luminance: inside
        // the Note 5f discontinuity window no nominal peak reaches the
        // requested value exactly, and the table has to agree with the
        // formulas rendered below it
        this.extendedPeakValue.textContent = `${data.peak109.toFixed(0)} cd/m²`;

        // Render the formulas from the solved configuration; the
        // Extended Range Formula doubles as the forward check of the
        // solved Lw against the requested extended peak
        this.formulaRenderer.updateFormulas(
            data.lw,
            surroundLuminance,
            blackLevel,
            blackLevelEnabled,
            data
        );
    }

    /**
     * Update extended range values based on current settings
     * 
     * @param {number} peakLuminance - Current peak luminance
     * @param {number} surroundLuminance - Current surround luminance
     * @param {number} blackLevel - Current black level
     * @param {boolean} blackLevelEnabled - Whether black level is enabled
     */
    updateExtendedRangeValues(peakLuminance, surroundLuminance, blackLevel = 0, blackLevelEnabled = false) {
        const effectiveBlackLevel = blackLevelEnabled ? blackLevel : 0;
    
        // Calculate extended range values using the simplified method
        const extendedRangeData = this.calculator.computeExtendedRangeValues(peakLuminance, surroundLuminance, effectiveBlackLevel);
    
        // Update the main table values
        this.uiController.displayDerivedValues(
            extendedRangeData.gamma,
            extendedRangeData.referenceWhite,
            extendedRangeData.blackLevelLift
        );
    
        // Update the displayed values
        this.extendedPeakValue.textContent = `${extendedRangeData.peak109.toFixed(0)} cd/m²`;
    
        // Update formulas with extended range data
        this.formulaRenderer.updateFormulas(
            peakLuminance,
            surroundLuminance,
            blackLevel,
            blackLevelEnabled,
            extendedRangeData
        );
    }
}

//=============================================================================
// MAIN APPLICATION CLASS
//=============================================================================

/**
 * Main application class that coordinates all components
 */
class HLGAdaptationCalculator {
    constructor() {
        // Initialize core components
        this.calculator = new HLGCalculator();
        this.uiController = new UIController(this.calculator);
        this.formulaRenderer = new FormulaRenderer(this.calculator);
        this.presetController = new PresetController(this.uiController);
        this.blackLevelController = new BlackLevelController(this.uiController);

        // Initialize optional components
        this.extendedRangeController = new ExtendedRangeController(
            this.calculator,
            this.uiController,
            this.formulaRenderer
        );

        this.sliderResetController = new SliderResetController(
            this.calculator,
            this.uiController,
            this.extendedRangeController
        );

        this.graphPanelController = new GraphPanelController(
            this.calculator,
            this.uiController,
            this.extendedRangeController
        );

        // In Extended Range mode the ExtendedRangeController drives the
        // formula rendering and (in extended-peak mode) the derived
        // table cells, so those components defer to it
        this.formulaRenderer.extendedRangeController = this.extendedRangeController;
        this.uiController.extendedRangeController = this.extendedRangeController;
    }

    /**
     * Initialize the application
     */
    initialize() {
        this.ensureKatexLoaded(() => {
            // Initialize all components
            this.presetController.initialize();
            this.blackLevelController.initialize();
            this.formulaRenderer.initialize();
            this.uiController.setupEventListeners();
            this.extendedRangeController.initialize();
            this.sliderResetController.initialize();
            this.graphPanelController.initialize();

            // Initial update
            this.uiController.updateDisplay(true);
            this.uiController.emitUpdateEvent();
        });
    }

    /**
     * Check if KaTeX is loaded, with retry
     * 
     * @param {Function} callback - Function to call when KaTeX is loaded
     * @param {number} attempts - Number of attempts so far
     */
    ensureKatexLoaded(callback, attempts = 0) {
        const maxAttempts = 20; // 20 * 50ms = 1 second max wait

        if (window.katex) {
            console.log("KaTeX loaded successfully");
            callback();
        } else if (attempts < maxAttempts) {
            console.log(`Waiting for KaTeX to load... (${attempts + 1}/${maxAttempts})`);
            setTimeout(() => this.ensureKatexLoaded(callback, attempts + 1), 50);
        } else {
            console.warn("KaTeX failed to load after multiple attempts. Formulas may not render correctly.");
            callback();
        }
    }
}

// Initialize the application when the document is loaded
document.addEventListener("DOMContentLoaded", () => {
    const app = new HLGAdaptationCalculator();
    app.initialize();
});
