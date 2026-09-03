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
 * The "Block Diagram and Graphs" panel of the HLG Display Adaptation
 * Calculator: the four graph views, the display controls they repeat
 * and the panel that owns them, drawn with D3.
 *
 * - CameraOetfGraphController: HLG Reference OETF graph view
 * - EotfGraphController: Live display EOTF curve graph view
 * - SurroundGammaGraphController: System Gamma vs surround graph view
 * - PeakGammaGraphController: System Gamma vs peak luminance graph view
 * - GraphDisplaySettings: The calculator's display controls repeated in a view
 * - GraphPanelController: d3 loading and switching between graph views
 *
 * Loaded before hlg_calculator.js, whose HLGAdaptationCalculator
 * constructs the GraphPanelController with the calculator and UI
 * controller the views read from.
 */

//=============================================================================
// GRAPH LABEL SYMBOLS
//=============================================================================

/**
 * Set an SVG text element's content with its quantity symbols
 * italic, as BT.2100-3 Table 5 typesets them and as <var> does in
 * the panel's HTML. SVG has no <var>, so a symbol is marked with
 * asterisks in the string and drawn as its own tspan:
 * "Signal level *E\u2032* (%)". The content is rebuilt on each
 * call, the labels carrying a live value being re-set as it moves.
 */
/**
 * The probe's mark on a curve: a reticle, where the levels a chart
 * states are plain dots. Its centre is a dot of the markers' own
 * size, so it reads as the same kind of point and the arms are all
 * that say the reader moves this one; they begin clear of the dot's
 * white ring. They run with the axes, so the arm along the view's
 * crosshair extends that line and the other reads as a tick across
 * it. Every view's probe uses it, so the reader learns one shape.
 */
const PROBE_ARMS = "M-12,0L-7.5,0M7.5,0L12,0M0,-12L0,-7.5M0,7.5L0,12";

/**
 * Add a probe reticle, drawn about its own origin so a translate
 * places it and a display style hides it
 */
function appendProbeMark(parent, role) {
    const mark = parent.append("g").attr("class", `hover-mark ${role}`);
    mark.append("path").attr("d", PROBE_ARMS);
    mark.append("circle").attr("r", 4.5);
    return mark;
}

function drawLabel(text, string) {
    text.text(null);
    string.split("*").forEach((part, i) => {
        if (part) {
            text.append("tspan")
                .attr("class", i % 2 ? "sym" : null)
                .text(part);
        }
    });
}

//=============================================================================
// GRAPH RESPONSIVE LAYOUT
//=============================================================================

/**
 * Rebuild a graph controller's chart at its container's actual pixel
 * width whenever that width changes — including each time the
 * container becomes visible: a hidden view measures 0 and is skipped,
 * and its text measures as nothing, so anything drawn while hidden is
 * rebuilt on showing. Charts are drawn in real pixels rather than scaled
 * from a fixed-size drawing, so text and touch targets keep their
 * designed size on phones and desktops alike.
 *
 * The controller must provide container, builtWidth, applyLayout(),
 * buildFrame() and redraw().
 */
function makeGraphResponsive(controller) {
    const rebuild = () => {
        const width = controller.container.clientWidth;
        if (width < 10) {
            controller.builtWidth = 0;
            return;
        }
        if (Math.abs(width - controller.builtWidth) < 2) return;
        controller.builtWidth = width;
        controller.applyLayout(width);
        controller.container.innerHTML = "";
        controller.buildFrame();
        controller.redraw();
    };

    if (window.ResizeObserver) {
        new ResizeObserver(rebuild).observe(controller.container);
    } else {
        window.addEventListener("resize", rebuild);
        rebuild();
    }
}

/**
 * Wire a graph probe's slider and value box to a setter, as the
 * calculator's own sliders are wired: either control moves the
 * probe and keeps the other in step. A box being typed into is
 * left alone (setValue gets fromInput = true, so a half-typed
 * "0." survives) and re-formats on leaving.
 *
 * A range input stops only on multiples of its step, so its
 * maximum is unreachable by dragging unless it is one — and both
 * probes end at the 109% signal level, (1019-64)/(940-64) =
 * 1.090182648 (code 1019, Table 9), which is not. The last step
 * therefore reports the exact maximum, so the probe lands where
 * the Extended Range Formula does. An untouched box is never read
 * back either: it shows a rounded value ("109.02"), which would
 * quantize the probe to the box's own precision.
 *
 * sliderToValue converts the slider's units to the probe's, for
 * the System Gamma views, whose sliders move in log10(cd/m2) —
 * a linear slider would spend most of its travel in the top
 * decade — while their boxes stay in cd/m2.
 */
function wireGraphProbeControls(slider, input, setValue, sliderToValue = v => v) {
    const snapToMax = value => {
        const max = parseFloat(slider.max);
        const step = parseFloat(slider.step);
        return (Number.isFinite(max) && Number.isFinite(step) && max - value < step) ? max : value;
    };
    let edited = false;
    slider.addEventListener("input", () => {
        edited = false;
        setValue(sliderToValue(snapToMax(parseFloat(slider.value))), false);
    });
    input.addEventListener("input", () => {
        edited = true;
        const typed = parseFloat(input.value);
        if (Number.isFinite(typed)) setValue(typed, true);
    });
    input.addEventListener("blur", () => {
        if (!edited) return;
        edited = false;
        const typed = parseFloat(input.value);
        setValue(Number.isFinite(typed) ? typed : 0, false);
    });
}

//=============================================================================
// CAMERA OETF GRAPH CONTROLLER
//=============================================================================

/**
 * HLG Reference OETF view — the camera curve the system is named
 * after (BT.2100-3 Table 5): the square-root segment, approximating
 * SDR cameras, joined at the E' = 0.5 breakpoint to the logarithmic
 * segment that extends the highlight range, with the super-white
 * continuation above 100% dashed (Note 5h). Behind it runs an SDR
 * camera at the same sensitivity, the BT.709-6 OETF (see
 * sdrSignalAt), which clips at 1.15x Reference White where HLG
 * reaches 3.77x — the comparison BT.2390-12 section 6.1 draws in
 * its Figure 18.
 *
 * Table 5 fixes E'(E) but not the scene light behind E: Note 5b
 * leaves that mapping free, so the scene axis is a choice of
 * exposure, set above the chart as the other views set their
 * conditions. The same grey card reads 38%, 21.2% or 67.2%
 * depending on the mapping.
 *
 * No calculator setting reaches the curve; the exposure does, so
 * setExposure rebuilds it. The probe converts scene light to
 * signal level and 10-bit code value (Table 9), camera direction.
 */
class CameraOetfGraphController {
    constructor(calculator, uiController, extendedRangeController) {
        this.calculator = calculator;

        this.container = document.getElementById("cameraOetfGraph");
        this.readout = document.getElementById("cameraOetfReadout");

        // Persistent probe: crosshair + readout position in scene
        // light, driven by the slider and value box under the
        // chart. Defaults to the white card (scene 1.0).
        this.PROBE_DEFAULT = 1.0;
        this.probeSlider = document.getElementById("cameraOetfProbe");
        this.probeInput = document.getElementById("cameraOetfProbeInput");
        this.probeResetButton = document.getElementById("cameraOetfProbeResetBtn");
        this.probeScene = this.PROBE_DEFAULT;

        // The exposure setting over the chart: slider, value box
        // and preset row (see buildExposurePresets)
        this.exposureSlider = document.getElementById("cameraOetfExposure");
        this.exposureInput = document.getElementById("cameraOetfExposureInput");
        this.resetButton = document.getElementById("cameraOetfResetBtn");
        this.presetRow = document.getElementById("cameraOetfPresets");
        this.presetButtons = [];

        // The SDR trace's legend entry, which the exposure renames
        this.sdrDetail = document.getElementById("cameraOetfSdrDetail");

        // The anchor line under the chart, carrying what the
        // exposure moves at each of the chart's markers
        this.anchorLine = document.getElementById("cameraOetfAnchors");

        // Chart geometry; set from the container's real width by
        // applyLayout() before every (re)build
        this.width = 720;
        this.height = 420;
        this.margin = { top: 30, right: 20, bottom: 52, left: 56 };
        this.narrow = false;
        this.builtWidth = 0;

        // The camera exposure: the Table 5 E value a 100% white
        // card produces. Note 5b leaves the mapping free, so it is
        // this view's one control rather than a constant. Starts at
        // the BT.2408-9 line-up, HDR Reference White at 75%.
        this.EXPOSURE_DEFAULT = this.calculator.hlgInverseOetf(0.75);
        this.exposure = this.EXPOSURE_DEFAULT;

        // The line-ups Note 5b's freedom has been exercised into,
        // in the order the slider runs: WHP 309 / ARIB STD-B67
        // v1.0's anchoring, which is also the bottom of its range,
        // BT.2408-9's, and BT.2390-12 Figure 18's equal-sensitivity
        // setup. The range runs on above them to an exposure no
        // document proposes as a line-up, which reports as custom.
        this.PRESETS = [
            [1 / 12, "ARIB STD-B67 v1.0"],
            [this.EXPOSURE_DEFAULT, "ITU-R BT.2408-9"],
            [1 / 3, "ITU-R BT.2390-12 Fig 18"]
        ];

        // What the exposure may be set to, as the signal level of
        // the 18% grey card: from the WHP 309 / ARIB STD-B67 v1.0
        // anchoring at one end to feeding scene light to Table 5
        // unconverted at the other, which are the two ways the same
        // card is read at 21% and at 67%.
        this.GREY_LEVEL_MIN = this.calculator.hlgOetf(0.18 / 12) * 100;
        this.GREY_LEVEL_MAX = this.calculator.hlgOetf(0.18) * 100;

        // The axis steps between these, taking the smallest that
        // still holds the 109% ceiling. That ceiling moves twelve
        // times over the exposure range — 6.19x the white card at
        // the BT.2408-9 line-up, 19.68x at ARIB STD-B67 v1.0's —
        // and no one scale serves both: 6.4 loses the ceiling
        // below about a 37% grey card, 20 leaves the default view's
        // landmarks in the left third. Tracking the curve instead
        // would hold its shape fixed and move only the tick
        // numbers, so the slider would appear to do nothing.
        // 6.4 is the default exposure's scale.
        this.SCENE_GREY = 0.18;
        this.SCENE_MIN = 0;
        this.X_SCALES = [2, 4, 6.4, 10, 20];
        this.X_TICKS = {
            2: [0, 0.5, 1, 1.5, 2],
            4: [0, 1, 2, 3, 4],
            6.4: [0, 1, 2, 3, 4, 5, 6],
            10: [0, 2, 4, 6, 8, 10],
            20: [0, 4, 8, 12, 16, 20]
        };

        this.SAMPLES = 160;

        this.recomputeLandmarks();
    }

    /**
     * Recompute what the current exposure puts where.
     *
     * Scene light is held in multiples of a 100% reflectance
     * white card, so
     * the cards never move — their 2.47 stops apart is a
     * reflectance ratio, not a camera setting. What a *signal*
     * level defines moves instead: the breakpoint, the 100% and
     * 109% ceilings, and where the SDR camera clips.
     */
    recomputeLandmarks() {
        this.SCENE_BREAK = (1 / 12) / this.exposure;
        this.SCENE_100 = 1 / this.exposure;
        this.SCENE_109 = this.calculator.hlgInverseOetf(
            this.calculator.SIGNAL_LEVEL_109) / this.exposure;

        // The SDR comparison camera's exposure, matched to HLG's on
        // the grey card, and where it clips (see sdrSignalAt)
        this.SDR_EXPOSURE = this.sdrInverseOetf(this.greyLevel()) / this.SCENE_GREY;
        this.SCENE_SDR_CLIP = 1 / this.SDR_EXPOSURE;

        // The smallest scale that still holds the 109% ceiling
        this.SCENE_MAX = this.X_SCALES.find(m => m >= this.SCENE_109)
            || this.X_SCALES[this.X_SCALES.length - 1];

        // The probe's ceiling: the exact 109% scene level, which
        // wireGraphProbeControls' last step lands on (it is not a
        // multiple of the controls' step)
        this.SCENE_PROBE_MAX = Math.min(this.SCENE_109, this.SCENE_MAX);
    }

    /**
     * The signal levels the exposure puts the two cards at. The
     * white card is HDR Reference White when it lands on 75%
     * (BT.2408-9 section 2.1)
     */
    greyLevel(exposure = this.exposure) {
        return this.calculator.hlgOetf(this.SCENE_GREY * exposure);
    }

    whiteLevel() {
        return this.calculator.hlgOetf(this.exposure);
    }

    /**
     * Set the exposure from the grey card's signal level, the
     * quantity BT.2390-12 section 6.1 configures a camera by.
     * Redraws: the curve itself depends on it
     */
    setExposure(greyPercent, fromInput = false) {
        const level = Math.max(this.GREY_LEVEL_MIN,
            Math.min(this.GREY_LEVEL_MAX, greyPercent)) / 100;
        this.exposure = this.calculator.hlgInverseOetf(level) / this.SCENE_GREY;
        this.recomputeLandmarks();
        // The box is left alone while it is the one being typed in
        if (!fromInput) this.syncExposureInput();
        this.rebuildAtCurrentWidth();
    }

    /**
     * Set a preset's exposure itself rather than the grey level the
     * slider works in, so updateLegend matches it outright instead
     * of within the tolerance a rounded level would need
     */
    setExposureExact(exposure) {
        this.exposure = exposure;
        this.recomputeLandmarks();
        this.syncExposureInput();
        this.rebuildAtCurrentWidth();
    }

    syncExposureInput() {
        this.exposureInput.value = (this.greyLevel() * 100).toFixed(2);
    }

    /**
     * Fill the preset row. Each button carries the level its
     * exposure puts the grey card at, which is what makes the four
     * readable against one another and against the box beside them.
     */
    buildExposurePresets() {
        this.presetRow.innerHTML = "";
        this.presetButtons = this.PRESETS.map(([exposure, name]) => {
            const button = document.createElement("button");
            const level = document.createElement("span");
            level.className = "preset-level";
            level.textContent = `${(this.greyLevel(exposure) * 100).toFixed(1)} %`;
            button.type = "button";
            // The issuer is a span of its own, for the phone row to
            // drop: a BT number cites without it
            const issuer = name.startsWith("ITU-R ") ? "ITU-R " : "";
            if (issuer) {
                const span = document.createElement("span");
                span.className = "preset-issuer";
                span.textContent = issuer;
                button.append(span);
            }
            button.append(name.slice(issuer.length), level);
            button.onclick = () => this.setExposureExact(exposure);
            this.presetRow.append(button);
            return button;
        });
    }

    /**
     * Rebuild at the width already measured. makeGraphResponsive
     * rebuilds only when the width changes, which is the right
     * test for a resize and the wrong one for an exposure change
     */
    rebuildAtCurrentWidth() {
        if (!this.builtWidth) return;
        this.applyLayout(this.builtWidth);
        this.container.innerHTML = "";
        this.buildFrame();
        this.redraw();
    }

    /**
     * Initialize the graph (called by GraphPanelController once
     * d3 is known to be available). The chart itself is built by
     * makeGraphResponsive() as soon as the container has a width.
     */
    initialize() {
        wireGraphProbeControls(this.probeSlider, this.probeInput,
            (scene, fromInput) => this.setProbe(scene, fromInput));
        this.probeResetButton.onclick = () => this.setProbe(this.PROBE_DEFAULT);

        wireGraphProbeControls(this.exposureSlider, this.exposureInput,
            (greyPercent, fromInput) => this.setExposure(greyPercent, fromInput));
        this.resetButton.onclick = () => this.setExposureExact(this.EXPOSURE_DEFAULT);
        this.buildExposurePresets();
        makeGraphResponsive(this);
    }

    /**
     * Compute geometry for the given container width (narrow layouts
     * get a fixed height and tighter margins)
     */
    applyLayout(width) {
        this.narrow = width < 520;
        this.width = width;
        this.height = this.narrow ? 340 : Math.min(420, Math.round(width * 0.55));
        this.margin = { top: 30, right: this.narrow ? 14 : 20, bottom: 52, left: 56 };
    }

    /**
     * Signal level E' for a scene light value in white-card units
     */
    signalAt(scene) {
        return this.calculator.hlgOetf(scene * this.exposure);
    }

    /**
     * The SDR comparison camera: the BT.709-6 OETF (item 1.2),
     * V = 1.099 L^0.45 − 0.099 above L = 0.018 and 4.5 L below,
     * which BT.2020 shares and BT.2390-12 Figure 18 draws beside
     * HLG. The linear toe is the difference section 6.1 notes
     * below 50% signal, HLG using the square root throughout.
     *
     * "Same sensitivity" is the Report's own definition: both
     * cameras give the same signal on the 18% grey card. The HLG
     * exposure puts the card at a level, the SDR exposure is
     * whatever puts BT.709's V there, and L = 1 then clips it.
     *
     * The Report's printed levels — "SDR (89%)" and "SDR (100%)"
     * in Figure 19, the "factor of 3" — follow from idealizing
     * this camera as a pure square root, which puts its clip at a
     * third of the 100% HLG scene light; the BT.709 curve clips
     * earlier.
     */
    sdrOetf(L) {
        return L < 0.018 ? 4.5 * L : 1.099 * Math.pow(L, 0.45) - 0.099;
    }

    sdrInverseOetf(V) {
        return V < 4.5 * 0.018 ? V / 4.5 : Math.pow((V + 0.099) / 1.099, 1 / 0.45);
    }

    sdrSignalAt(scene) {
        return this.sdrOetf(Math.min(1, scene * this.SDR_EXPOSURE));
    }

    /**
     * Build the chart. Everything is static: the OETF does not
     * depend on the calculator's inputs.
     */
    buildFrame() {
        const { width, height, margin } = this;

        this.svg = d3.select(this.container).append("svg")
            .attr("width", width)
            .attr("height", height)
            .classed("narrow", this.narrow)
            .attr("role", "img")
            .attr("aria-label", "The HLG Reference OETF per BT.2100-3 Table 5: signal level against linear scene light, with the square-root and logarithmic segments marked");

        this.x = d3.scaleLinear()
            .domain([this.SCENE_MIN, this.SCENE_MAX])
            .range([margin.left, width - margin.right]);
        this.y = d3.scaleLinear()
            .domain([0, 1.15])
            .range([height - margin.bottom, margin.top]);

        // Horizontal grid at the signal-level ticks, to the 109%
        // super-white limit. The scale stays in Table 5's E' units,
        // which the curve is computed in; only the labels are the
        // percent the rest of the view speaks
        const yTicks = [0, 0.25, 0.5, 0.75, 1.0, 1.09];
        this.svg.append("g").attr("class", "grid")
            .selectAll("line")
            .data(yTicks)
            .join("line")
            .attr("x1", margin.left).attr("x2", width - margin.right)
            .attr("y1", d => this.y(d)).attr("y2", d => this.y(d));

        // Axes
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(0,${height - margin.bottom})`)
            .call(d3.axisBottom(this.x)
                .tickValues(this.X_TICKS[this.SCENE_MAX])
                .tickFormat(v => String(v))
                .tickSizeOuter(0));
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(${margin.left},0)`)
            .call(d3.axisLeft(this.y)
                .tickValues(yTicks)
                .tickFormat(v => (v * 100).toFixed(0))
                .tickSizeOuter(0));

        // Axis titles
        this.svg.append("text").attr("class", "axis-title")
            .attr("x", (margin.left + width - margin.right) / 2)
            .attr("y", height - 10)
            .attr("text-anchor", "middle")
            .text("Linear scene light (1.0 = 100 % reflectance white card)");
        this.svg.append("text").attr("class", "axis-title")
            .attr("transform", `translate(16, ${(margin.top + height - margin.bottom) / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .call(t => drawLabel(t, "Signal level *E\u2032* (%)"));

        // Vertical markers at the signal-level anchors: the 18%
        // grey card's nominal 38% (BT.2408-9 Table 1, scene 0.18),
        // the white card (75% at that line-up), 100% and 109%
        //
        // The exposure moves the two kinds in opposite ways: a card
        // keeps its place on the axis and its label moves, while a
        // coding level is fixed and the exposure moves where in the
        // scene it falls — off the axis entirely, at a low enough
        // exposure, and then it is not drawn.
        const level = v => `${(v * 100).toFixed(1)} %`;
        // A coding level's label carries the scene light it takes to
        // reach it, which is the headroom a camera needs above the
        // white card; a card's label cannot, its place being fixed
        const ceiling = (name, scene) =>
            this.narrow ? name : `${name} (${scene.toFixed(2)}\u00d7)`;
        const markers = [
            [this.SCENE_GREY, level(this.greyLevel())],
            [1.0, level(this.whiteLevel())],
            [this.SCENE_100, ceiling("100 %", this.SCENE_100)],
            [this.SCENE_109, ceiling("109 %", this.SCENE_109)]
        ].filter(([scene], i, all) =>
            scene >= this.SCENE_MIN && scene <= this.SCENE_MAX
            // At the top of the exposure range the card and the
            // 100% signal are one point, and one line cannot carry
            // two labels; the coding level keeps its own, being the
            // landmark the card has landed on. The tolerance is
            // float noise, not a spacing rule: markers merely close
            // together are what the staggered rows below are for.
            && !all.some(([other], j) => j > i && Math.abs(other - scene) < 1e-9));
        // Stagger the labels by parity. Stops put 100% and 109%
        // within a stop of each other, and the exposure can bring
        // any two markers together — at the top of its range the
        // white card and the 100% signal are the same point — so
        // the rows alternate at every width, not just on a phone
        markers.forEach(([scene, label], i) => {
            // A ceiling sits close to the right-hand edge by
            // construction — the axis is scaled to hold it — so a
            // label centred on one would overhang the chart. Those
            // anchor to the edge instead; the half-width is estimated
            // from the label, the text not being measurable until the
            // view is shown, and the anchor set as a style, the
            // class carrying one of its own.
            const px = this.x(scene);
            const edge = width - margin.right;
            const overhang = px + label.length * 3.4 > edge;
            this.svg.append("line")
                .attr("class", "ref-line")
                .attr("x1", px).attr("x2", px)
                .attr("y1", margin.top).attr("y2", height - margin.bottom);
            this.svg.append("text")
                .attr("class", "ref-line-label keep-narrow")
                .attr("x", overhang ? edge : px)
                .style("text-anchor", overhang ? "end" : null)
                .attr("y", margin.top - (i % 2 ? 20 : 8))
                .text(label);
        });

        // The curve, in its named segments: square-root ("gamma") up
        // to the E' = 0.5 breakpoint, logarithmic to the 100% signal,
        // dashed continuation through the super-whites
        // The middle of the logarithmic segment as drawn, which the
        // exposure moves and the axis' edge can cut short
        const logMid = (this.SCENE_BREAK
            + Math.min(this.SCENE_100, this.SCENE_MAX)) / 2;

        const line = d3.line()
            .x(d => this.x(d[0]))
            .y(d => this.y(d[1]));
        const sample = (s0, s1) => d3.range(this.SAMPLES + 1).map(i => {
            const scene = s0 + (s1 - s0) * i / this.SAMPLES;
            return [scene, this.signalAt(scene)];
        });

        // What the exposure has carried off the right-hand edge is
        // not drawn, so every segment stops at the axis
        const onAxis = v => Math.min(v, this.SCENE_MAX);

        // SDR comparison trace (BT.709-6 item 1.2; see
        // sdrSignalAt). Drawn first, so the HLG curve sits over it
        // where the two nearly coincide, below the breakpoint
        this.svg.append("path").attr("class", "curve-sdr")
            .attr("d", line(sample(this.SCENE_MIN, this.SCENE_MAX)
                .map(([scene]) => [scene, this.sdrSignalAt(scene)])));

        this.svg.append("path").attr("class", "curve-oetf-gamma")
            .attr("d", line(sample(this.SCENE_MIN, onAxis(this.SCENE_BREAK))));
        if (this.SCENE_BREAK < this.SCENE_MAX) {
            this.svg.append("path").attr("class", "curve-oetf-log")
                .attr("d", line(sample(this.SCENE_BREAK, onAxis(this.SCENE_100))));
        }
        if (this.SCENE_100 < this.SCENE_MAX) {
            this.svg.append("path").attr("class", "curve-oetf-super")
                .attr("d", line(sample(this.SCENE_100, onAxis(this.SCENE_109))));
        }

        // Breakpoint dot (scene 1/12 in the OETF's normalization)
        if (this.SCENE_BREAK < this.SCENE_MAX) {
            this.svg.append("circle")
                .attr("class", "marker-dot")
                .attr("r", 4.5)
                .attr("cx", this.x(this.SCENE_BREAK))
                .attr("cy", this.y(0.5));
        }


        // Direct labels on leader lines. A label naming something
        // the exposure moves is positioned against that landmark,
        // or it would drift off what it names; the rest are
        // fractions of the axis, which changes scale under them.
        const across = f => this.SCENE_MAX * f;
        const labels = [
            {
                text: this.narrow ? "Square-root segment" : "Square-root (gamma) segment \u2014 approximates SDR cameras",
                tx: across(0.1172), ty: 0.30,
                cx: this.SCENE_BREAK / 2, cy: this.signalAt(this.SCENE_BREAK / 2)
            },
            {
                text: this.narrow ? "Log segment" : "Logarithmic segment \u2014 extends highlight range",
                tx: across(0.3984), ty: 0.68,
                cx: logMid, cy: this.signalAt(logMid)
            },
            {
                // Named in both units: the axis reads percent, but
                // the join is where E' = 0.5, the value-continuity
                // condition c is solved from. In ARIB STD-B67 v1.0
                // that point is reference white (r = 0.5, code 502),
                // and at the BT.2408-9 exposure this chart carries
                // both it and the white card, so the dot is named
                text: this.narrow
                    ? "Breakpoint"
                    : "50 % breakpoint (*E\u2032* = 0.5) \u2014 reference white in ARIB STD-B67 v1.0",
                tx: across(0.1172), ty: 0.52,
                cx: this.SCENE_BREAK + across(0.00625), cy: 0.5
            },
            {
                // Sits directly above the SDR clip plateau, so it
                // needs no leader line
                text: this.narrow
                    ? `SDR clips at ${this.SCENE_SDR_CLIP.toFixed(2)}\u00d7`
                    : `SDR \u2014 same exposure, clips at ${this.SCENE_SDR_CLIP.toFixed(2)}\u00d7`,
                tx: this.SCENE_SDR_CLIP + across(0.03125), ty: 1.045, noLeader: true
            },
            {
                // Sits under the dashed continuation above the 100%
                // signal (BT.2100-3 Note 5h), the only curve up there
                text: "Super-whites",
                tx: this.SCENE_100 + across(0.0547), ty: 0.93, noLeader: true
            }
        ].filter(label => label.tx < this.SCENE_MAX
            && (label.cx === undefined || label.cx < this.SCENE_MAX));

        labels.forEach(label => {
            if (!label.noLeader) {
                this.svg.append("line")
                    .attr("class", "segment-line")
                    .attr("x1", this.x(label.tx) - 4).attr("y1", this.y(label.ty))
                    .attr("x2", this.x(label.cx)).attr("y2", this.y(label.cy));
            }
            this.svg.append("text")
                .attr("class", "segment-label")
                .attr("x", this.x(label.tx)).attr("y", this.y(label.ty) + 4)
                .call(t => drawLabel(t, label.text));
        });

        // Probe layer: crosshair and a dot on the curve (the values
        // go to the HTML readout line below the chart)
        this.hoverGroup = this.svg.append("g")
            .attr("class", "hover")
            .style("display", "none");
        this.hoverLine = this.hoverGroup.append("line")
            .attr("class", "hover-line")
            .attr("y1", margin.top)
            .attr("y2", height - margin.bottom);
        this.hoverDot = appendProbeMark(this.hoverGroup, "dot-adapted");

        this.exposureSlider.min = this.GREY_LEVEL_MIN.toFixed(2);
        this.exposureSlider.max = this.GREY_LEVEL_MAX.toFixed(2);
        this.exposureInput.min = this.exposureSlider.min;
        this.exposureInput.max = this.exposureSlider.max;
        this.syncExposureSlider();
        this.updateLegend();

        // The probe controls span the axis' scene light range
        this.probeSlider.min = this.SCENE_MIN;
        this.probeSlider.max = this.SCENE_PROBE_MAX;
        this.probeInput.min = this.SCENE_MIN;
        this.probeInput.max = this.SCENE_PROBE_MAX;
    }

    /**
     * Keep the slider on the exposure's level, to its own step,
     * when the level was set some other way
     */
    syncExposureSlider() {
        const level = (this.greyLevel() * 100).toFixed(1);
        if (Math.abs(parseFloat(this.exposureSlider.value) - level) > 0.05) {
            this.exposureSlider.value = level;
        }
    }

    /**
     * Whether the exposure is set at a given one of the PRESETS
     */
    setAt(target) {
        return Math.abs(this.exposure - target) < target * 0.002;
    }

    /**
     * Report the exposure: which of the PRESETS is in force, marked
     * on the row itself; the level the SDR trace puts diffuse
     * white at; and the anchor line under
     * the chart, with what the exposure moves at each marker.
     */
    updateLegend() {
        const setAt = target => this.setAt(target);
        const named = this.PRESETS.find(([g]) => setAt(g));
        this.presetButtons.forEach((button, i) => {
            const active = this.PRESETS[i] === named;
            button.classList.toggle("active", active);
            button.setAttribute("aria-pressed", String(active));
        });

        // The SDR trace's entry carries the level it puts the white
        // card at, as BT.2390-12 Figure 19 names its traces, and
        // where it clips, which is what the trace is drawn to show.
        // Once the exposure carries the card past the clip it has
        // no level to report.
        const clip = `clips at ${this.SCENE_SDR_CLIP.toFixed(2)}\u00d7 white`;
        this.sdrDetail.textContent = this.SCENE_SDR_CLIP >= 1
            ? ` (ITU-R BT.709-6 OETF; white card at ${(this.sdrSignalAt(1) * 100).toFixed(1)} %, ${clip})`
            : ` (ITU-R BT.709-6 OETF; ${clip}, below the white card)`;

        // The chart's four markers, each with the value the exposure
        // moves: a card keeps its scene light and its signal level
        // moves, a coding level keeps its signal level and the scene
        // light it takes moves. Every value names its axis, the view
        // being about percentages that do not say which axis they
        // are of. The grey card is the input, being what the
        // exposure is set as; the rest follow from it.
        //
        // The card is what the exposure moves; the level it lands
        // on is HDR Reference White only where BT.2408-9 section
        // 2.1 puts it, so only there does the level take that name.
        //
        // A coding level's scene light is quoted as BT.2390-12
        // section 6.1 quotes it, "scene luminance of 375% diffuse
        // white" and "about 620% if super-whites are used", with
        // the multiplier the top axis shows alongside. It is the
        // unit a camera's dynamic range is set in — percent of the
        // white it is lined up to, the quantity NHK's 8K HLG
        // cameras are set by (Funatsu et al., IS&T Electronic
        // Imaging 2017, p. 153: sensor saturation as a percentage
        // of reference white, a 1200 % setting filling the signal
        // to 100 % and about 2000 % to 109 %, which the ARIB
        // line-up reproduces here) — and what the sensor has to
        // reach for the signal to use that level at all.
        const level = value => `${(value * 100).toFixed(1)} % signal`;
        const ofWhite = scene =>
            `${(scene * 100).toFixed(0)} % of diffuse white (${scene.toFixed(2)}\u00d7)`;
        const anchors = [
            ["18 % grey card", level(this.greyLevel()), "input"],
            [setAt(this.EXPOSURE_DEFAULT) ? "HDR Reference White" : "white card",
                level(this.whiteLevel()), "calculated"],
            ["100 % signal", ofWhite(this.SCENE_100), "calculated"],
            ["109 % signal", ofWhite(this.SCENE_109), "calculated"]
        ];
        const ROLE_LABELS = { input: "Input", calculated: "Calculated" };
        this.anchorLine.innerHTML =
            anchors.map(([name, value, role]) =>
                `<span class="readout-item readout-${role}">` +
                `<span class="anchor-level">${name}</span>` +
                `<span class="anchor-value">${value}</span>` +
                `<span class="role-pill role-${role}">${ROLE_LABELS[role]}</span>` +
                `</span>`
            ).join(" ");
    }

    /**
     * The curve is static; redraw re-renders the probe (the chart
     * is rebuilt on resize)
     */
    redraw() {
        if (!this.svg) return;
        this.setProbe(this.probeScene);
    }

    /**
     * Move the persistent probe to a scene light value. fromInput
     * leaves the value box's text alone while it is being typed in.
     */
    setProbe(scene, fromInput = false) {
        this.probeScene = Math.max(this.SCENE_MIN, Math.min(this.SCENE_PROBE_MAX, scene));
        this.probeSlider.value = this.probeScene;
        if (!fromInput) this.probeInput.value = this.probeScene.toFixed(2);
        this.renderProbe();
    }

    /**
     * Render the crosshair, curve dot and readout at the probe
     * position: scene light \u2192 signal level \u2192 10-bit code value
     */
    renderProbe() {
        if (!this.svg) return;

        const scene = this.probeScene;
        const signal = this.signalAt(scene);
        const code = Math.round(64 + signal * 876);

        this.hoverGroup.style("display", null);
        this.hoverLine.attr("x1", this.x(scene)).attr("x2", this.x(scene));
        this.hoverDot.attr("transform", `translate(${this.x(scene)},${this.y(signal)})`);

        this.readout.innerHTML =
            `<span class="readout-item">Scene light ${scene.toFixed(2)}\u00d7` +
            ` (${(scene * 100).toFixed(0)} % of white card)</span> \u2192 ` +
            `<span class="readout-item"><var>E\u2032</var> = ${(signal * 100).toFixed(2)}%</span> \u2192 ` +
            `<span class="readout-item">CV ${code}</span>`;
    }
}

//=============================================================================
// EOTF GRAPH CONTROLLER
//=============================================================================

/**
 * Live plot of the HLG Reference EOTF (BT.2100-3 Table 5): signal
 * level E' against the display luminance it produces, for the
 * current configuration — system gamma per Note 5f and BT.2390-12
 * Section 6.2, black level lift when enabled — against the fixed
 * reference display (Lw = 1000 cd/m2, 5 cd/m2 surround, γ = 1.2).
 * The display adaptation is the gap between the two curves.
 *
 * Signal on the vertical axis and light on the horizontal, the
 * orientation a waveform monitor puts the signal in. The
 * luminance axis is logarithmic and fixed: fitted to the data it
 * would swallow the motion the settings cause. The right-hand
 * axis pairs the signal with its 10-bit code values (Table 9:
 * code = 876·E' + 64), as a scope's own scale does. The probe is
 * driven by a slider rather than hover, so touch and keyboard
 * reach it.
 *
 * Display-side throughout: the EOTF takes a signal and a display,
 * and the camera exposure is no term of it. So nothing here moves
 * with the OETF view's exposure — every level marked is fixed, and
 * what the display settings change is the luminance each reaches.
 */
class EotfGraphController {
    constructor(calculator, uiController, extendedRangeController) {
        this.calculator = calculator;
        this.uiController = uiController;
        this.extendedRangeController = extendedRangeController;

        this.container = document.getElementById("eotfGraph");

        this.readout = document.getElementById("eotfReadout");

        // The current display's anchor luminances, listed under the
        // chart; rebuilt with the markers on every redraw
        this.anchorLine = document.getElementById("eotfAnchors");

        // The system gamma the display settings above the chart
        // produce, stated under it
        this.gammaLine = document.getElementById("eotfGamma");

        // The range mode's pill, which survives the narrow layout's
        // hiding of the legend detail
        this.rangeMode = document.getElementById("eotfRangeMode");

        // BT.2408-9 section 2.2 Table 1's two nominal levels, the
        // reference values this view is drawn against. Both are fixed:
        // Table 1's levels "do not change" with the display, and none
        // of them is the camera's to move. The grey card's follows
        // from Reference White by geometry — 18 % of it in scene
        // light — which is where BT.2408-9's rounded 38 % comes from;
        // the exact level is kept, so the luminance it reaches on the
        // reference display comes out at Table 1's own 26 cd/m².
        this.WHITE_LEVEL = 75;
        this.GREY_LEVEL = calculator.hlgOetf(
            0.18 * calculator.hlgInverseOetf(0.75)) * 100;

        // Persistent probe: crosshair + readout position as a signal
        // fraction, driven by the slider and value box under the
        // chart. Defaults to HDR Reference White (75%).
        this.PROBE_DEFAULT = this.WHITE_LEVEL / 100;
        this.probeSlider = document.getElementById("eotfProbe");
        this.probeInput = document.getElementById("eotfProbeInput");
        this.probeResetButton = document.getElementById("eotfProbeResetBtn");
        this.probeSignal = this.PROBE_DEFAULT;

        // Chart geometry; set from the container's real width by
        // applyLayout() before every (re)build
        this.width = 720;
        this.height = 420;
        this.margin = { top: 28, right: 62, bottom: 56, left: 52 };
        this.narrow = false;
        this.builtWidth = 0;

        // Fixed luminance domain (see class comment). The floor clips
        // signal levels whose luminance falls below 0.001 cd/m²; the
        // ceiling covers the extended-peak maximum of 25000 cd/m².
        this.L_MIN = 0.001;
        this.L_MAX = 30000;
        this.L_TICKS = [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000];
        // Narrow layouts label every other decade; the grid keeps them all
        this.L_TICKS_NARROW = [0.001, 0.1, 10, 1000];

        this.SAMPLES = 240;

        // Parameters of the currently drawn curve, kept for the hover readout
        this.current = null;
    }

    /**
     * Initialize the graph (called by GraphPanelController once
     * d3 is known to be available). The chart itself is built by
     * makeGraphResponsive() as soon as the container has a width.
     */
    initialize() {
        this.setupEventListeners();
        makeGraphResponsive(this);
    }

    /**
     * Compute geometry for the given container width (narrow layouts
     * get a fixed height and tighter margins)
     */
    applyLayout(width) {
        this.narrow = width < 520;
        this.width = width;
        // Taller than the other views by the 1.17 the R 103 envelope
        // costs: the nominal video range is 100 of the 116.8 points of
        // signal on the axis, and this leaves its curve the room the
        // other views give theirs
        this.height = this.narrow ? 396 : Math.min(490, Math.round(width * 0.64));
        this.margin = { top: 28, right: this.narrow ? 46 : 62, bottom: 56, left: this.narrow ? 42 : 52 };
    }

    /**
     * Build the static SVG frame: scales, grid, axes, curve
     * paths, marker and hover layers. The luminance axis and grid
     * never change; redraw() only updates the signal and code axes,
     * the band backdrop, the curve paths and the markers.
     */
    buildFrame() {
        const { width, height, margin } = this;

        this.svg = d3.select(this.container).append("svg")
            .attr("width", width)
            .attr("height", height)
            .classed("narrow", this.narrow)
            .attr("role", "img")
            .attr("aria-label", "The HLG Reference EOTF per BT.2100-3 Table 5: signal level versus the display luminance it produces, for the current settings and the reference display");

        // Luminance runs along the bottom, signal up the side (see
        // the class comment); the signal domain is set in redraw()
        this.x = d3.scaleLog()
            .domain([this.L_MIN, this.L_MAX])
            .range([margin.left, width - margin.right]);
        this.y = d3.scaleLinear()
            .range([height - margin.bottom, margin.top]);

        // Clip the curves to the plot area (the log floor cuts them off)
        this.svg.append("defs").append("clipPath")
            .attr("id", "eotfClip")
            .append("rect")
            .attr("x", margin.left).attr("y", margin.top)
            .attr("width", width - margin.left - margin.right)
            .attr("height", height - margin.top - margin.bottom);

        // Signal-range band backdrop (filled in redraw(): the bands
        // depend on the range mode), behind the grid
        this.bandGroup = this.svg.append("g");

        // Vertical grid at the luminance decades
        this.svg.append("g").attr("class", "grid")
            .selectAll("line")
            .data(this.L_TICKS)
            .join("line")
            .attr("y1", margin.top).attr("y2", height - margin.bottom)
            .attr("x1", d => this.x(d)).attr("x2", d => this.x(d));

        // The bands' own line and names go in front of the grid,
        // which would otherwise rule through them
        this.bandOverlayGroup = this.svg.append("g");

        // Axes: the luminance axis is fixed, the two signal axes are
        // set in redraw()
        const luminanceTickFormat = v =>
            v >= 1000 ? v.toLocaleString("en-US").replace(/,/g, " ") : String(v);
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(0,${height - margin.bottom})`)
            .call(d3.axisBottom(this.x)
                .tickValues(this.narrow ? this.L_TICKS_NARROW : this.L_TICKS)
                .tickFormat(luminanceTickFormat)
                .tickSizeOuter(0));
        this.yAxisGroup = this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(${margin.left},0)`);

        // Second signal axis on the right: the same levels as 10-bit
        // narrow-range code values (BT.2100-3 Table 9: code =
        // 876·E' + 64), the pairing a waveform monitor's own scale
        // carries; tick values set in redraw()
        this.codeAxisGroup = this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(${width - margin.right},0)`);

        // Axis titles
        this.svg.append("text").attr("class", "axis-title")
            .attr("x", (margin.left + width - margin.right) / 2)
            .attr("y", height - 10)
            .attr("text-anchor", "middle")
            .text("Display luminance (cd/m², log scale)");
        this.svg.append("text").attr("class", "axis-title")
            .attr("transform", `translate(14, ${(margin.top + height - margin.bottom) / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .call(t => drawLabel(t, "Signal level *E′* (%)"));
        this.svg.append("text").attr("class", "axis-title")
            .attr("transform", `translate(${width - 12}, ${(margin.top + height - margin.bottom) / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .text("10-bit code value (CV)");

        // Curves (reference behind the adapted curve)
        const curves = this.svg.append("g").attr("clip-path", "url(#eotfClip)");
        this.referencePath = curves.append("path").attr("class", "curve-reference");
        this.adaptedPath = curves.append("path").attr("class", "curve-adapted");

        this.markerGroup = this.svg.append("g");

        // Probe layer: crosshair and a dot on each curve (the values
        // go to the HTML readout line above the chart)
        this.hoverGroup = this.svg.append("g")
            .attr("class", "hover")
            .style("display", "none");
        // The probe is a signal level, so its crosshair runs across
        // the chart at that level
        this.hoverLine = this.hoverGroup.append("line")
            .attr("class", "hover-line")
            .attr("x1", margin.left)
            .attr("x2", width - margin.right);
        this.hoverDotReference = appendProbeMark(this.hoverGroup, "dot-reference");
        this.hoverDotAdapted = appendProbeMark(this.hoverGroup, "dot-adapted");

    }

    /**
     * Set up redraw triggers and the probe slider
     */
    setupEventListeners() {
        document.addEventListener('hlg-values-updated', () => this.redraw());

        wireGraphProbeControls(this.probeSlider, this.probeInput,
            (percent, fromInput) => this.setProbe(percent / 100, fromInput));
        this.probeResetButton.onclick = () => this.setProbe(this.PROBE_DEFAULT);

        // The range-mode buttons re-present without emitting
        // hlg-values-updated. These listeners are registered after the
        // ExtendedRangeController's, so its mode flags are already
        // updated when the graph redraws.
        document.getElementById("nominalRangeButton").addEventListener("click", () => this.redraw());
        document.getElementById("extendedRangeButton").addEventListener("click", () => this.redraw());

    }

    /**
     * Display luminance (cd/m²) of the currently drawn adapted
     * configuration at signal level E' (achromatic, R'=G'=B'=E')
     */
    adaptedLuminanceAt(signal) {
        const { lw, gamma, beta, extendedRange } = this.current;
        if (!extendedRange) signal = Math.min(signal, 1.0); // super-whites clip at 100%
        const [luminance] = this.calculator.hlgEotf(signal, signal, signal, gamma, lw, 0, beta);
        return luminance;
    }

    /**
     * Display luminance (cd/m²) of the reference display
     * (Lw = 1000 cd/m², 5 cd/m² surround, γ = 1.2) at signal level E'
     */
    referenceLuminanceAt(signal) {
        if (!this.current.extendedRange) signal = Math.min(signal, 1.0); // super-whites clip at 100%
        const gamma = this.calculator.calculateSystemGamma(this.calculator.L_REF, this.calculator.LS_REF);
        const [luminance] = this.calculator.hlgEotf(signal, signal, signal, gamma, this.calculator.L_REF, 0, null);
        return luminance;
    }

    /**
     * Redraw the curves and markers from the current UI state
     */
    redraw() {
        if (!this.svg) return;

        const { peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled } = this.uiController.getCurrentValues();
        const effectiveBlackLevel = blackLevelEnabled ? blackLevel : 0;
        const extendedRange = this.extendedRangeController.extendedRangeActive;

        // In extended-peak mode the slider value is the extended peak;
        // the curve must be drawn from the solved nominal Lw (the same
        // computation the ExtendedRangeController presents in the table)
        const data = (extendedRange && this.extendedRangeController.sliderSetsExtendedPeak)
            ? this.calculator.solveNominalFromExtendedPeak(peakLuminance, surroundLuminance, effectiveBlackLevel)
            : this.calculator.computeExtendedRangeValues(peakLuminance, surroundLuminance, effectiveBlackLevel);

        // The signal axis always spans the full coding range to 109%:
        // Nominal Range mode does not shrink the axis, it clips —
        // signals above 100% display at the nominal peak (the flat
        // top of the curve) and the super-white band greys out
        const maxSignal = this.calculator.SIGNAL_LEVEL_109;
        this.current = {
            lw: data.lw,
            gamma: data.gamma,
            beta: data.blackLevelLift > 0 ? data.blackLevelLift : null,
            maxSignal,
            extendedRange
        };

        // The current display's own conditions are the controls above
        // the chart and the values under it — the two settings, the
        // system gamma they produce (Note 5f, BT.2390-12 section 6.2)
        // and the anchor luminances — so the legend entry carries only
        // the range mode, as a pill on the curve's name. The reference
        // display's conditions appear nowhere else, so its entry keeps
        // them in the markup
        this.rangeMode.textContent = extendedRange ? "Extended Range" : "Nominal Range";
        this.rangeMode.className =
            `legend-mode legend-mode-${extendedRange ? "extended" : "nominal"}`;

        // In extended-peak mode γ follows the solved nominal peak, not
        // the peak the slider above sets
        this.gammaLine.innerHTML =
            `<span class="readout-item readout-calculated">System Gamma <var>γ</var> = ${data.gamma.toFixed(3)}` +
            `<span class="role-pill role-calculated">Calculated</span></span>`;

        // The signal axis carries the whole 10-bit code range, so
        // that R 103 Figure 1's bands can be drawn in full: the code
        // values SDI reserves bound it top and bottom (E' = 0 is CV
        // 64 and 100% is CV 940, so E' = (CV - 64) / 876)
        const atCode = cv => (cv - 64) / 876 * 100;
        this.y.domain([atCode(0), atCode(1023)]);
        this.yAxisGroup.call(d3.axisLeft(this.y)
            .tickValues([atCode(4), 0, 25, 50, 75, 100, 109])
            .tickFormat(d => d === atCode(4) ? "-6.84" : d)
            .tickSizeOuter(0));

        const codeTicks = [4, 64, 721, 940, 1019].map(atCode);
        this.codeAxisGroup.call(d3.axisRight(this.y)
            .tickValues(codeTicks)
            .tickFormat(pct => Math.round(64 + 876 * pct / 100))
            .tickSizeOuter(0));

        // Band backdrop, in the colours and terms of EBU R 103
        // Figure 1: nominal video range (CV 64-940), the headroom
        // either side (CV 941-1019, CV 4-63), and SDI's Time
        // Reference Signal values (CV 1020-1023, CV 0-3). Headroom
        // stays yellow in both range modes, as the Figure has it;
        // only the curve and the band's label change
        this.bandGroup.selectAll("*").remove();
        const bandX = this.margin.left;
        const bandW = this.width - this.margin.left - this.margin.right;
        const band = (c0, c1, cls) => this.bandGroup.append("rect")
            .attr("class", cls)
            .attr("x", bandX).attr("y", this.y(atCode(c1)))
            .attr("width", bandW).attr("height", this.y(atCode(c0)) - this.y(atCode(c1)));
        band(0, 4, "zone-novideo");
        band(4, 64, "zone-headroom");
        band(64, 940, "zone-nominal");
        band(940, 1019, "zone-headroom");
        band(1019, 1023, "zone-novideo");

        this.bandOverlayGroup.selectAll("*").remove();
        // R 103's two preferred limits, each tagged as the Figure
        // tags them
        [105, -5].forEach(pct => {
            this.bandOverlayGroup.append("line")
                .attr("class", "ref-line")
                .attr("x1", bandX).attr("x2", bandX + bandW)
                .attr("y1", this.y(pct)).attr("y2", this.y(pct));
            // Each label sits on the side of its limit facing the
            // nominal video range: the outer side of the -5% line is
            // only a couple of code values deep
            this.bandOverlayGroup.append("text")
                .attr("class", "pct-label")
                .attr("x", bandX + bandW / 2).attr("y", this.y(pct) + (pct > 0 ? 13 : -5))
                .attr("text-anchor", "middle")
                .text(this.narrow ? `${pct} %` : `EBU R 103, ${pct} %`);
        });
        // R 103 Figure 1 names its red bands; here the axis is to
        // scale and they are four code values tall, with no room —
        // the code axis ending at 1019 and 4 says the same. Each
        // headroom band is named in the corner the curve leaves
        // empty, on its own R 103 limit's line: centred, the name
        // would sit a couple of pixels off that dashed line
        this.bandOverlayGroup.append("text")
            .attr("class", "band-label")
            .attr("x", bandX + 6).attr("y", this.y(105) + 13)
            .text(extendedRange ? "headroom" : "headroom, clipped");
        this.bandOverlayGroup.append("text")
            .attr("class", "band-label")
            .attr("x", bandX + 6).attr("y", this.y(-5) - 5)
            .text("headroom");
        // Turned along the band it names, and centred on it: the
        // nominal range is an extent on the signal axis, where the
        // headroom tags above and below name bands too thin to hold
        // anything but a corner
        this.bandOverlayGroup.append("text")
            .attr("class", "band-label")
            .attr("transform", `translate(${bandX + bandW - 8}, ${(this.y(0) + this.y(100)) / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .text("Nominal video range");

        // BT.2408-9's two reference levels, marked across the plot.
        // They are a matched pair from one table, so they are named
        // the same way; what the display settings change is the
        // luminance each reaches, not the level. The labels sit on the
        // left, the side of the plot the curve has left by these
        // levels.
        const levelLabels = [
            [this.WHITE_LEVEL, "HDR Reference White", "HDR Ref White"],
            [this.GREY_LEVEL, "18 % grey card"]
        ].map(([level, name, shortName = name]) => {
            this.bandOverlayGroup.append("line")
                .attr("class", "level-line")
                .attr("x1", bandX).attr("x2", bandX + bandW)
                .attr("y1", this.y(level)).attr("y2", this.y(level));
            const text = this.bandOverlayGroup.append("text")
                .attr("class", "pct-label")
                .attr("x", bandX + 6).attr("y", this.y(level) - 5)
                .text(this.narrow
                    ? shortName
                    : `${name}, ${level.toFixed(level < 50 ? 1 : 0)} % (ITU-R BT.2408-9)`);
            return { level, text, shortName };
        });

        // Each curve starts exactly where it crosses the log floor,
        // solved by bisection, or uniform sampling leaves a visible
        // gap at the steep bottom end. The EOTF is monotonic; black
        // lift can put E' = 0 above the floor, and it starts at 0
        const sample = luminanceAt => {
            let s0 = 0;
            if (luminanceAt(0) < this.L_MIN) {
                let lo = 0, hi = maxSignal;
                for (let i = 0; i < 40; i++) {
                    const mid = (lo + hi) / 2;
                    if (luminanceAt(mid) < this.L_MIN) lo = mid; else hi = mid;
                }
                s0 = hi;
            }
            return d3.range(this.SAMPLES + 1).map(i => {
                const signal = s0 + (maxSignal - s0) * i / this.SAMPLES;
                return [signal, luminanceAt(signal)];
            });
        };
        const line = d3.line()
            .defined(d => d[1] >= this.L_MIN)
            .x(d => this.x(d[1]))
            .y(d => this.y(d[0] * 100));
        this.adaptedPath.attr("d", line(sample(s => this.adaptedLuminanceAt(s))));
        this.referencePath.attr("d", line(sample(s => this.referenceLuminanceAt(s))));

        // Markers on the adapted curve: HDR Reference White, nominal
        // peak, and — in Extended Range mode, where super-whites are
        // displayed rather than clipped — the extended peak. Each
        // carries the results table's colour role: the level the
        // peak slider sets is the input, the rest calculated, so the
        // 100% and 109% roles swap with the peak mode
        const setsExtendedPeak = extendedRange && this.extendedRangeController.sliderSetsExtendedPeak;
        const markers = [
            { signal: this.GREY_LEVEL / 100, label: "18 % grey card", role: "calculated" },
            { signal: this.WHITE_LEVEL / 100, label: "75%", role: "calculated" },
            { signal: 1.0, label: "100%", role: setsExtendedPeak ? "calculated" : "input" }
        ];
        if (extendedRange) {
            markers.push({
                signal: this.calculator.SIGNAL_LEVEL_109,
                label: "109%",
                role: setsExtendedPeak ? "input" : "calculated"
            });
        }

        this.markerGroup.selectAll("*").remove();
        const markerLabels = [];
        markers.forEach(marker => {
            const luminance = this.adaptedLuminanceAt(marker.signal);
            if (luminance < this.L_MIN) return;

            const cx = this.x(luminance);
            const group = this.markerGroup.append("g")
                .attr("transform", `translate(${cx},${this.y(marker.signal * 100)})`);
            group.append("circle").attr("class", "marker-dot").attr("r", 4.5);
            // The luminance the level reaches, which is what the
            // display settings change and what no axis gives directly
            // — a log decade hides where 203 falls between 100 and
            // 1 000. The level itself is on the y axis, and named
            // across the plot for the two the standards name.
            //
            // The label sits left of its dot, on the empty side of a
            // curve that climbs to the right. The 109% label clears
            // the 100% one on every layout: the nine signal points
            // between the dots are 11 px on the narrow plot, more
            // than the label's height, and the top margin holds it.
            const text = group.append("text").attr("class", "marker-label label-left")
                .attr("x", -9)
                .attr("y", -8)
                .text(`${this.formatLuminance(luminance)} cd/m²`);
            markerLabels.push({ level: marker.signal * 100, cx, text });
        });

        // Use short name to avoid overlapping with luminance label
        levelLabels.forEach(({ level, text, shortName }) => {
            const marker = markerLabels.find(m => Math.abs(m.level - level) < 1e-9);
            if (!marker) return;
            const nameEnd = bandX + 6 + text.node().getComputedTextLength();
            const labelStart = marker.cx - 9 - marker.text.node().getComputedTextLength();
            if (nameEnd + 8 > labelStart) text.text(shortName);
        });

        // The same anchors listed under the chart with their
        // luminances, which the curve labels have no room for and
        // the probe reaches one at a time. Each carries the results
        // table's role colour and pill, in the table's order, so a
        // level reads the same in both places. No curve swatch: blue
        // here is the input role, not the curve.
        //
        // Level, value, pill is the order the phone layout stacks;
        // the stylesheet puts the pill back between the pair on
        // wider screens, so neither layout needs a wrapper
        const ROLE_LABELS = { input: "Input", calculated: "Calculated" };
        // The grey card leads them, the level being the lowest. It is
        // named rather than given as a level: the card is what it
        // identifies, and its line on the chart carries the percentage
        const anchors = markers;
        this.anchorLine.innerHTML =
            anchors.map(anchor =>
                `<span class="readout-item readout-${anchor.role}">` +
                `<span class="anchor-level">${anchor.label}</span>` +
                `<span class="anchor-value">${this.formatLuminance(this.adaptedLuminanceAt(anchor.signal))} cd/m²</span>` +
                `<span class="role-pill role-${anchor.role}">${ROLE_LABELS[anchor.role]}</span>` +
                `</span>`
            ).join(" ");

        // Re-render the probe against the new curve (setProbe
        // re-clamps and keeps both controls in step). The slider
        // carries the exact endpoint, 109.018265%, so its last step
        // lands on it; the value box carries it as displayed,
        // 109.02%, so typing what the box shows is in range
        this.probeSlider.max = maxSignal * 100;
        this.probeInput.max = (maxSignal * 100).toFixed(2);
        this.setProbe(this.probeSignal);
    }

    /**
     * Format a luminance value for the hover readout
     */
    formatLuminance(luminance) {
        if (luminance >= 100) return luminance.toFixed(0);
        if (luminance >= 1) return luminance.toFixed(1);
        if (luminance >= 0.01) return luminance.toFixed(3);
        return luminance.toFixed(4);
    }

    /**
     * Move the persistent probe to a signal fraction. fromInput
     * leaves the value box's text alone while it is being typed in.
     */
    setProbe(signal, fromInput = false) {
        if (!this.current) return;
        this.probeSignal = Math.max(0, Math.min(this.current.maxSignal, signal));
        this.probeSlider.value = this.probeSignal * 100;
        // Two decimals: one 10-bit code value is 1/876 = 0.114% of
        // signal, so a coarser box could not tell the 109% endpoint
        // (109.02%, code 1019) from the 109.0% one step below it
        if (!fromInput) this.probeInput.value = (this.probeSignal * 100).toFixed(2);
        this.renderProbe();
    }

    /**
     * Render the crosshair, curve dots and readout at the probe
     * position. The probe persists and is re-rendered on every
     * redraw, so it tracks the live settings.
     */
    renderProbe() {
        if (!this.current) return;

        const signal = this.probeSignal;
        const py = this.y(signal * 100);

        const adapted = this.adaptedLuminanceAt(signal);
        const reference = this.referenceLuminanceAt(signal);

        this.hoverGroup.style("display", null);
        this.hoverLine.attr("y1", py).attr("y2", py);

        // dy lifts a label clear of the crosshair or drops it below;
        // the x is clamped so a reading near either end of the axis
        // stays inside the plot
        const placeMark = (mark, luminance) => {
            if (luminance >= this.L_MIN) {
                mark.style("display", null)
                    .attr("transform", `translate(${this.x(luminance)},${py})`);
            } else {
                mark.style("display", "none");
            }
        };
        placeMark(this.hoverDotAdapted, adapted);
        placeMark(this.hoverDotReference, reference);

        // In Nominal Range mode the signal clips, not the luminance:
        // a probe in the super-white band is clamped to 1.0 before
        // the EOTF, so both displays read their nominal peak. The
        // note sits on the signal, never on a luminance it would
        // contradict
        const clipped = !this.current.extendedRange && signal > 1.0;
        // Last of the three, so it reads against the signal as a whole:
        // between the two it would bind to the code value instead, and
        // CV 1019 is this probe's own 109.02%, not the 100% clipped to
        // (that is CV 940, on the chart's top axis). "at" rather than
        // "to" for the same reason — it sets up no equivalence
        const clipNote = clipped
            ? `<span class="readout-item">(clipped at 100%)</span> `
            : "";
        this.readout.innerHTML =
            `<span class="readout-item"><var>E′</var> = ${(signal * 100).toFixed(2)}%</span> ` +
            `<span class="readout-item">CV ${Math.round(64 + signal * 876)}</span> ` +
            clipNote +
            `<span class="readout-item"><span class="readout-swatch swatch-adapted"></span>Current display ${this.formatLuminance(adapted)} cd/m²</span> ` +
            `<span class="readout-item"><span class="readout-swatch swatch-reference"></span>Reference display ${this.formatLuminance(reference)} cd/m²</span>`;
    }
}

//=============================================================================
// SURROUND GAMMA GRAPH CONTROLLER
//=============================================================================

/**
 * System Gamma vs surround luminance view — the view of ITU-R
 * BT.2390-12 Figure 21 ("Graph of system gamma vs. ambient lighting
 * for a number of different screen luminances").
 *
 * Drawn from calculateSystemGamma, not the Figure's lines of best
 * fit: BT.2390-12 Section 6.2 gives two surround models and this
 * tool implements the multiplicative one, so the curves differ
 * slightly from the Figure's subtractive fit.
 *
 * Grey context curves hold the nominal peak at each preset value
 * (400-10000 cd/m²), so every preset button has its curve; fixed,
 * they are drawn once in the frame. The blue curve follows the
 * current display's nominal Lw, solved in extended-peak mode. The
 * surround axis extends the Figure's 5-500 cd/m² down to 0.05.
 * γ does not depend on black level directly; Black Level Lift
 * reaches this view only in extended-peak mode, through the
 * nominal Lw the solver returns.
 */
class SurroundGammaGraphController {
    constructor(calculator, uiController, extendedRangeController) {
        this.calculator = calculator;
        this.uiController = uiController;
        this.extendedRangeController = extendedRangeController;

        this.container = document.getElementById("surroundGammaGraph");

        this.readout = document.getElementById("surroundGammaReadout");
        this.legendDetail = document.getElementById("surroundCurrentLegendDetail");

        // Persistent probe: crosshair + readout position in surround
        // luminance, driven by the slider and value box under the
        // chart. Defaults to the 5 cd/m² reference surround
        // (BT.2100-3 Table 3), the level the dashed line marks.
        this.PROBE_DEFAULT = calculator.LS_REF;
        this.probeSlider = document.getElementById("surroundGammaProbe");
        this.probeInput = document.getElementById("surroundGammaProbeInput");
        this.probeResetButton = document.getElementById("surroundGammaProbeResetBtn");
        this.probeLs = this.PROBE_DEFAULT;

        // Chart geometry; set from the container's real width by
        // applyLayout() before every (re)build. The wide right
        // margin holds the family end labels.
        this.width = 720;
        this.height = 420;
        this.margin = { top: 28, right: 64, bottom: 52, left: 70 };
        this.narrow = false;
        this.builtWidth = 0;

        // Surround axis spans the slider range; the gamma axis is
        // fixed so the operating point visibly moves (headroom above
        // 2.0 covers the extended-peak mode maximum)
        this.LS_MIN = 0.05;
        this.LS_MAX = 500;
        this.X_TICKS = [0.05, 0.5, 5, 50, 500];
        this.GAMMA_MIN = 0.6;
        this.GAMMA_MAX = 2.1;
        this.Y_TICKS = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];

        // Fixed nominal peaks of the context curves: the peak
        // luminance preset values (kept in sync with
        // PresetController.PEAK_LUMINANCE_PRESETS), so every preset
        // button has its curve on the chart
        this.FAMILY_LW = [400, 500, 600, 1000, 2000, 3000, 4000, 5000, 10000];

        this.SAMPLES = 160;

        // Parameters of the currently drawn curve, kept for the hover readout
        this.current = null;
    }

    /**
     * Initialize the graph (called by GraphPanelController once
     * d3 is known to be available). The chart itself is built by
     * makeGraphResponsive() as soon as the container has a width.
     */
    initialize() {
        this.setupEventListeners();
        makeGraphResponsive(this);
    }

    /**
     * Compute geometry for the given container width (narrow layouts
     * get a fixed height and a tighter label margin)
     */
    applyLayout(width) {
        this.narrow = width < 520;
        this.width = width;
        this.height = this.narrow ? 340 : Math.min(420, Math.round(width * 0.55));
        this.margin = { top: 28, right: this.narrow ? 50 : 64, bottom: 52, left: 70 };
    }

    /**
     * Build the static SVG frame. Everything except the
     * current-display curve, its operating point and the hover layer
     * is fixed and drawn once here, including the family curves.
     */
    buildFrame() {
        const { width, height, margin } = this;

        this.svg = d3.select(this.container).append("svg")
            .attr("width", width)
            .attr("height", height)
            .classed("narrow", this.narrow)
            .attr("role", "img")
            .attr("aria-label", "System gamma versus surround luminance for the current display and for fixed nominal peak luminances");

        this.x = d3.scaleLog()
            .domain([this.LS_MIN, this.LS_MAX])
            .range([margin.left, width - margin.right]);
        this.y = d3.scaleLinear()
            .domain([this.GAMMA_MIN, this.GAMMA_MAX])
            .range([height - margin.bottom, margin.top]);

        this.svg.append("defs").append("clipPath")
            .attr("id", "surroundGammaClip")
            .append("rect")
            .attr("x", margin.left).attr("y", margin.top)
            .attr("width", width - margin.left - margin.right)
            .attr("height", height - margin.top - margin.bottom);

        // Horizontal grid at the gamma ticks
        this.svg.append("g").attr("class", "grid")
            .selectAll("line")
            .data(this.Y_TICKS)
            .join("line")
            .attr("x1", margin.left).attr("x2", width - margin.right)
            .attr("y1", d => this.y(d)).attr("y2", d => this.y(d));

        // Reference surround marker at 5 cd/m²
        this.svg.append("line")
            .attr("class", "ref-line")
            .attr("x1", this.x(5)).attr("x2", this.x(5))
            .attr("y1", margin.top).attr("y2", height - margin.bottom);
        this.svg.append("text")
            .attr("class", "ref-line-label")
            .attr("x", this.x(5)).attr("y", margin.top - 8)
            .text("reference surround");

        // Axes (both fixed)
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(0,${height - margin.bottom})`)
            .call(d3.axisBottom(this.x)
                .tickValues(this.X_TICKS)
                .tickFormat(v => String(v))
                .tickSizeOuter(0));
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(${margin.left},0)`)
            .call(d3.axisLeft(this.y)
                .tickValues(this.Y_TICKS)
                .tickFormat(v => v.toFixed(1))
                .tickSizeOuter(0));

        // Axis titles
        this.svg.append("text").attr("class", "axis-title")
            .attr("x", (margin.left + width - margin.right) / 2)
            .attr("y", height - 10)
            .attr("text-anchor", "middle")
            .text("Surround luminance (cd/m², log scale)");
        this.svg.append("text").attr("class", "axis-title")
            .attr("transform", `translate(16, ${(margin.top + height - margin.bottom) / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .text("System gamma");

        // Family curves for the Figure 21 nominal peaks, with the
        // peak value as an end label
        const curves = this.svg.append("g").attr("clip-path", "url(#surroundGammaClip)");
        this.FAMILY_LW.forEach(lw => {
            curves.append("path")
                .attr("class", "curve-family")
                .attr("d", this.lineGenerator()(this.sampleCurveAtPeak(lw)));
        });

        // End labels, nudged apart where neighbouring curves end too
        // close together (the 400/500/600 and 3000/4000/5000 preset
        // clusters end within a label height of each other)
        const MIN_LABEL_GAP = 12;
        const labels = this.FAMILY_LW.map(lw => ({
            lw,
            y: this.y(this.calculator.calculateSystemGamma(lw, this.LS_MAX)) + 4
        })).sort((a, b) => a.y - b.y);
        labels.forEach((label, i) => {
            if (i > 0 && label.y < labels[i - 1].y + MIN_LABEL_GAP) {
                label.y = labels[i - 1].y + MIN_LABEL_GAP;
            }
        });
        labels.forEach(label => {
            this.svg.append("text")
                .attr("class", "family-label")
                .attr("x", width - margin.right + 6)
                .attr("y", label.y)
                .text(label.lw.toLocaleString("en-US").replace(/,/g, " "));
        });

        // Header for the family end labels: without it the label
        // column reads as a second y-axis
        this.svg.append("text")
            .attr("class", "family-label")
            .attr("x", width - 4)
            .attr("y", margin.top - 8)
            .attr("text-anchor", "end")
            .text("nominal peak luminance (cd/m²)");

        // Current-display curve and operating point
        this.currentPath = curves.append("path").attr("class", "curve-adapted");
        this.operatingPoint = this.svg.append("g");
        this.operatingPoint.append("circle").attr("class", "marker-dot").attr("r", 4.5);
        this.operatingLabel = this.operatingPoint.append("text")
            .attr("class", "marker-label")
            .attr("y", -10);

        // Probe layer: crosshair and a dot on the current and the
        // 1000 cd/m² family curve (the values go to the HTML
        // readout line below the chart)
        this.hoverGroup = this.svg.append("g")
            .attr("class", "hover")
            .style("display", "none");
        this.hoverLine = this.hoverGroup.append("line")
            .attr("class", "hover-line")
            .attr("y1", margin.top)
            .attr("y2", height - margin.bottom);
        this.hoverDotReference = appendProbeMark(this.hoverGroup, "dot-reference");
        this.hoverDotCurrent = appendProbeMark(this.hoverGroup, "dot-adapted");

        // The probe controls span the axis
        this.probeSlider.min = Math.log10(this.LS_MIN);
        this.probeSlider.max = Math.log10(this.LS_MAX);
        this.probeInput.min = this.LS_MIN;
        this.probeInput.max = this.LS_MAX;
    }

    /**
     * Set up redraw triggers (same triggers as the EOTF view)
     * and the probe controls
     */
    setupEventListeners() {
        document.addEventListener('hlg-values-updated', () => this.redraw());
        document.getElementById("nominalRangeButton").addEventListener("click", () => this.redraw());
        document.getElementById("extendedRangeButton").addEventListener("click", () => this.redraw());

        // The slider moves in log10(cd/m²) so its travel matches the
        // axis; the value box stays in cd/m²
        wireGraphProbeControls(this.probeSlider, this.probeInput,
            (ls, fromInput) => this.setProbe(ls, fromInput),
            logLs => Math.pow(10, logLs));
        this.probeResetButton.onclick = () => this.setProbe(this.PROBE_DEFAULT);
    }

    /**
     * Move the persistent probe to a surround luminance. fromInput
     * leaves the value box's text alone while it is being typed in.
     */
    setProbe(ls, fromInput = false) {
        this.probeLs = Math.max(this.LS_MIN, Math.min(this.LS_MAX, ls));
        this.probeSlider.value = Math.log10(this.probeLs);
        if (!fromInput) this.probeInput.value = this.formatSurround(this.probeLs);
        this.renderProbe();
    }

    /**
     * d3 line generator for γ(Ls) samples
     */
    lineGenerator() {
        return d3.line()
            .x(d => this.x(d[0]))
            .y(d => this.y(d[1]));
    }

    /**
     * Sample γ(Ls) across the surround axis (geometrically, matching
     * the log axis) for a fixed nominal peak luminance
     */
    sampleCurveAtPeak(lw) {
        return d3.range(this.SAMPLES + 1).map(i => {
            const ls = this.LS_MIN * Math.pow(this.LS_MAX / this.LS_MIN, i / this.SAMPLES);
            return [ls, this.calculator.calculateSystemGamma(lw, ls)];
        });
    }

    /**
     * Redraw the current-display curve and operating point
     */
    redraw() {
        if (!this.svg) return;

        const { peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled } = this.uiController.getCurrentValues();
        const extendedRange = this.extendedRangeController;

        // In extended-peak mode the slider value is the extended
        // peak; the curve must be drawn from the solved nominal Lw
        let lw = peakLuminance;
        if (extendedRange.extendedRangeActive && extendedRange.sliderSetsExtendedPeak) {
            const effectiveBlackLevel = blackLevelEnabled ? blackLevel : 0;
            lw = this.calculator.solveNominalFromExtendedPeak(peakLuminance, surroundLuminance, effectiveBlackLevel).lw;
        }

        const gamma = this.calculator.calculateSystemGamma(lw, surroundLuminance);
        this.current = { lw, ls: surroundLuminance, gamma };

        this.currentPath.attr("d", this.lineGenerator()(this.sampleCurveAtPeak(lw)));

        // Operating point with the current γ as its label; nudge the
        // label inside the plot near the edges
        const cx = this.x(surroundLuminance);
        const xRange = this.x.range();
        this.operatingPoint.attr("transform", `translate(${cx},${this.y(gamma)})`);
        this.operatingLabel
            .attr("x", Math.min(Math.max(0, xRange[0] + 30 - cx), xRange[1] - 30 - cx))
            .call(t => drawLabel(t, `*γ* = ${gamma.toFixed(3)}`));

        // Nominal peak is the fixed parameter here — surround is the
        // axis, and the dot marks the display's own surround on it —
        // so the legend carries the peak the blue curve is drawn for
        this.legendDetail.innerHTML =
            ` (nominal peak <var>L<sub>W</sub></var> = ${this.formatPeak(lw)} cd/m²)`;

        // The probe persists across redraws, so it tracks the live
        // settings (setProbe re-clamps and keeps both controls in step)
        this.setProbe(this.probeLs);
    }

    /**
     * Format a surround luminance value for the hover readout
     */
    formatSurround(ls) {
        if (ls >= 10) return ls.toFixed(0);
        if (ls >= 1) return ls.toFixed(1);
        return ls.toFixed(2);
    }

    /**
     * Format a nominal peak luminance value for the hover readout
     */
    formatPeak(lw) {
        return Math.round(lw).toLocaleString("en-US").replace(/,/g, " ");
    }

    /**
     * Render the crosshair, curve dots and readout at the probe
     * position. The probe persists and is re-rendered on every
     * redraw, so it tracks the live settings.
     */
    renderProbe() {
        if (!this.current) return;

        const ls = this.probeLs;
        const px = this.x(ls);

        const gammaCurrent = this.calculator.calculateSystemGamma(this.current.lw, ls);
        const gammaReference = this.calculator.calculateSystemGamma(this.calculator.L_REF, ls);

        this.hoverGroup.style("display", null);
        this.hoverLine.attr("x1", px).attr("x2", px);
        this.hoverDotCurrent.attr("transform", `translate(${px},${this.y(gammaCurrent)})`);
        this.hoverDotReference.attr("transform", `translate(${px},${this.y(gammaReference)})`);

        this.readout.innerHTML =
            `<span class="readout-item">At ${this.formatSurround(ls)} cd/m²</span> ` +
            `<span class="readout-item"><span class="readout-swatch swatch-adapted"></span>Current (${this.formatPeak(this.current.lw)} cd/m²) <var>γ</var> = ${gammaCurrent.toFixed(3)}</span> ` +
            `<span class="readout-item"><span class="readout-swatch swatch-reference"></span>Reference (1&#8239;000 cd/m²) <var>γ</var> = ${gammaReference.toFixed(3)}</span>`;
    }
}

//=============================================================================
// PEAK GAMMA GRAPH CONTROLLER
//=============================================================================

/**
 * System Gamma vs nominal peak luminance view — the view of ITU-R
 * BT.2390-12 Figure 20 ("Gamma value to match images for different
 * screen peak brightness").
 *
 * Both Note 5f formulas run the full axis as grey context curves:
 * the basic one dashed, the extended κ one solid. The blue curve is
 * what calculateSystemGamma selects — basic within the 400-2000
 * cd/m² monitoring range, extended outside it — with extra samples
 * at the boundaries so the switchover steps stay crisp.
 *
 * The current surround factor applies to every curve, so the
 * operating dot stays on the blue one and the chart shifts with the
 * surround slider; at the 5 cd/m² reference surround the curves
 * match the Figure's. In extended-peak mode the dot sits at the
 * solved nominal Lw, pinned to the axis edge if it falls below it.
 */
class PeakGammaGraphController {
    constructor(calculator, uiController, extendedRangeController) {
        this.calculator = calculator;
        this.uiController = uiController;
        this.extendedRangeController = extendedRangeController;

        this.container = document.getElementById("peakGammaGraph");

        this.readout = document.getElementById("peakGammaReadout");
        this.legendDetail = document.getElementById("peakCurrentLegendDetail");

        // Persistent probe: crosshair + readout position in nominal
        // peak luminance, driven by the slider and value box under
        // the chart. Defaults to 1 000 cd/m², where BT.2100-3
        // Note 5f puts γ = 1.2.
        this.PROBE_DEFAULT = calculator.L_REF;
        this.probeSlider = document.getElementById("peakGammaProbe");
        this.probeInput = document.getElementById("peakGammaProbeInput");
        this.probeResetButton = document.getElementById("peakGammaProbeResetBtn");
        this.probeLw = this.PROBE_DEFAULT;

        // Chart geometry; set from the container's real width by
        // applyLayout() before every (re)build
        this.width = 720;
        this.height = 420;
        this.margin = { top: 28, right: 24, bottom: 52, left: 70 };
        this.narrow = false;
        this.builtWidth = 0;

        // Peak axis spans the nominal slider range; the ticks include
        // the Note 5f formula switchover values. The gamma axis
        // matches the System Gamma vs surround view.
        this.LW_MIN = 100;
        this.LW_MAX = 10000;
        this.X_TICKS = [100, 400, 1000, 2000, 10000];
        this.GAMMA_MIN = 0.6;
        this.GAMMA_MAX = 2.1;
        this.Y_TICKS = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];

        // Note 5f: basic formula within 400-2000 cd/m², extended outside
        this.MONITORING_RANGE = [400, 2000];

        this.SAMPLES = 200;

        // Parameters of the currently drawn curves, kept for the hover readout
        this.current = null;
    }

    /**
     * Initialize the graph (called by GraphPanelController once
     * d3 is known to be available). The chart itself is built by
     * makeGraphResponsive() as soon as the container has a width.
     */
    initialize() {
        this.setupEventListeners();
        makeGraphResponsive(this);
    }

    /**
     * Compute geometry for the given container width (narrow layouts
     * get a fixed height and tighter margins)
     */
    applyLayout(width) {
        this.narrow = width < 520;
        this.width = width;
        this.height = this.narrow ? 340 : Math.min(420, Math.round(width * 0.55));
        this.margin = { top: 28, right: this.narrow ? 18 : 24, bottom: 52, left: 70 };
    }

    /**
     * Build the static SVG frame: scales, grid, axes, the
     * production monitoring range markers, curve paths, the
     * operating point and the hover layer. The curves themselves
     * are drawn in redraw() — they shift with the surround factor.
     */
    buildFrame() {
        const { width, height, margin } = this;

        this.svg = d3.select(this.container).append("svg")
            .attr("width", width)
            .attr("height", height)
            .classed("narrow", this.narrow)
            .attr("role", "img")
            .attr("aria-label", "System gamma versus nominal peak luminance for both BT.2100-3 Note 5f formulas and the implemented formula switch");

        this.x = d3.scaleLog()
            .domain([this.LW_MIN, this.LW_MAX])
            .range([margin.left, width - margin.right]);
        this.y = d3.scaleLinear()
            .domain([this.GAMMA_MIN, this.GAMMA_MAX])
            .range([height - margin.bottom, margin.top]);

        this.svg.append("defs").append("clipPath")
            .attr("id", "peakGammaClip")
            .append("rect")
            .attr("x", margin.left).attr("y", margin.top)
            .attr("width", width - margin.left - margin.right)
            .attr("height", height - margin.top - margin.bottom);

        // Horizontal grid at the gamma ticks
        this.svg.append("g").attr("class", "grid")
            .selectAll("line")
            .data(this.Y_TICKS)
            .join("line")
            .attr("x1", margin.left).attr("x2", width - margin.right)
            .attr("y1", d => this.y(d)).attr("y2", d => this.y(d));

        // Production monitoring range markers (Note 5f switchovers)
        this.MONITORING_RANGE.forEach(lw => {
            this.svg.append("line")
                .attr("class", "ref-line")
                .attr("x1", this.x(lw)).attr("x2", this.x(lw))
                .attr("y1", margin.top).attr("y2", height - margin.bottom);
        });
        this.svg.append("text")
            .attr("class", "ref-line-label")
            .attr("x", this.x(Math.sqrt(this.MONITORING_RANGE[0] * this.MONITORING_RANGE[1])))
            .attr("y", margin.top - 8)
            .text("usual production monitoring range");

        // Axes (both fixed)
        const luminanceTickFormat = v =>
            v >= 1000 ? v.toLocaleString("en-US").replace(/,/g, " ") : String(v);
        // Narrow layouts keep only the decade ticks; the 400/2000
        // boundaries stay visible as the dashed range markers
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(0,${height - margin.bottom})`)
            .call(d3.axisBottom(this.x)
                .tickValues(this.narrow ? [100, 1000, 10000] : this.X_TICKS)
                .tickFormat(luminanceTickFormat)
                .tickSizeOuter(0));
        this.svg.append("g")
            .attr("class", "axis")
            .attr("transform", `translate(${margin.left},0)`)
            .call(d3.axisLeft(this.y)
                .tickValues(this.Y_TICKS)
                .tickFormat(v => v.toFixed(1))
                .tickSizeOuter(0));

        // Axis titles
        this.svg.append("text").attr("class", "axis-title")
            .attr("x", (margin.left + width - margin.right) / 2)
            .attr("y", height - 10)
            .attr("text-anchor", "middle")
            .text("Nominal peak luminance (cd/m², log scale)");
        this.svg.append("text").attr("class", "axis-title")
            .attr("transform", `translate(16, ${(margin.top + height - margin.bottom) / 2}) rotate(-90)`)
            .attr("text-anchor", "middle")
            .text("System gamma");

        // Curves: the two Note 5f formulas behind the implemented curve
        const curves = this.svg.append("g").attr("clip-path", "url(#peakGammaClip)");
        this.basicPath = curves.append("path").attr("class", "curve-family-dashed");
        this.extendedPath = curves.append("path").attr("class", "curve-family");
        this.implementedPath = curves.append("path").attr("class", "curve-adapted");

        // Operating point with the current γ as its label
        this.operatingPoint = this.svg.append("g");
        this.operatingPoint.append("circle").attr("class", "marker-dot").attr("r", 4.5);
        this.operatingLabel = this.operatingPoint.append("text")
            .attr("class", "marker-label")
            .attr("y", -10);

        // Probe layer: crosshair, a dot on the implemented curve and
        // one on the formula not in use at that peak (the values go to
        // the HTML readout line below the chart)
        this.hoverGroup = this.svg.append("g")
            .attr("class", "hover")
            .style("display", "none");
        this.hoverLine = this.hoverGroup.append("line")
            .attr("class", "hover-line")
            .attr("y1", margin.top)
            .attr("y2", height - margin.bottom);
        this.hoverDotOther = appendProbeMark(this.hoverGroup, "dot-reference");
        this.hoverDotImplemented = appendProbeMark(this.hoverGroup, "dot-adapted");

        // The probe controls span the axis
        this.probeSlider.min = Math.log10(this.LW_MIN);
        this.probeSlider.max = Math.log10(this.LW_MAX);
        this.probeInput.min = this.LW_MIN;
        this.probeInput.max = this.LW_MAX;
    }

    /**
     * Set up redraw triggers (same triggers as the other views)
     * and the probe controls
     */
    setupEventListeners() {
        document.addEventListener('hlg-values-updated', () => this.redraw());
        document.getElementById("nominalRangeButton").addEventListener("click", () => this.redraw());
        document.getElementById("extendedRangeButton").addEventListener("click", () => this.redraw());

        // The slider moves in log10(cd/m²) so its travel matches the
        // axis; the value box stays in cd/m²
        wireGraphProbeControls(this.probeSlider, this.probeInput,
            (lw, fromInput) => this.setProbe(lw, fromInput),
            logLw => Math.pow(10, logLw));
        this.probeResetButton.onclick = () => this.setProbe(this.PROBE_DEFAULT);
    }

    /**
     * Move the persistent probe to a nominal peak luminance.
     * fromInput leaves the value box's text alone while it is being
     * typed in.
     */
    setProbe(lw, fromInput = false) {
        this.probeLw = Math.max(this.LW_MIN, Math.min(this.LW_MAX, lw));
        this.probeSlider.value = Math.log10(this.probeLw);
        if (!fromInput) this.probeInput.value = Math.round(this.probeLw);
        this.renderProbe();
    }

    /**
     * d3 line generator for γ(Lw) samples
     */
    lineGenerator() {
        return d3.line()
            .x(d => this.x(d[0]))
            .y(d => this.y(d[1]));
    }

    /**
     * The surround factor of the multiplicative model, from the
     * calculator's constants (μ^log2(Ls/Ls_ref))
     */
    surroundFactor(ls) {
        const c = this.calculator;
        return Math.pow(c.MU, Math.log2(ls / c.LS_REF));
    }

    /**
     * The two Note 5f formulas drawn across the full axis for
     * comparison. These mirror the two branches of
     * HLGCalculator.calculateSystemGamma, which selects between
     * them; here each is needed on its own.
     */
    basicFormulaGamma(lw, ls) {
        const c = this.calculator;
        return (c.GAMMA_REF + 0.42 * Math.log10(lw / c.L_REF)) * this.surroundFactor(ls);
    }

    extendedRangeGamma(lw, ls) {
        const c = this.calculator;
        return c.GAMMA_REF * Math.pow(c.KAPPA, Math.log2(lw / c.L_REF)) * this.surroundFactor(ls);
    }

    /**
     * Sample γ(Lw) geometrically across the peak axis
     */
    sampleCurveFrom(gammaAt) {
        return d3.range(this.SAMPLES + 1).map(i => {
            const lw = this.LW_MIN * Math.pow(this.LW_MAX / this.LW_MIN, i / this.SAMPLES);
            return [lw, gammaAt(lw)];
        });
    }

    /**
     * Sample the implemented (piecewise) curve, adding sample pairs
     * at the switchover boundaries so the steps render vertically
     */
    sampleImplementedCurve(ls) {
        const lws = d3.range(this.SAMPLES + 1).map(i =>
            this.LW_MIN * Math.pow(this.LW_MAX / this.LW_MIN, i / this.SAMPLES));
        this.MONITORING_RANGE.forEach(boundary => {
            lws.push(boundary * 0.99999, boundary, boundary * 1.00001);
        });
        lws.sort((a, b) => a - b);
        return lws.map(lw => [lw, this.calculator.calculateSystemGamma(lw, ls)]);
    }

    /**
     * Redraw all curves and the operating point
     */
    redraw() {
        if (!this.svg) return;

        const { peakLuminance, surroundLuminance, blackLevel, blackLevelEnabled } = this.uiController.getCurrentValues();
        const extendedRange = this.extendedRangeController;

        // In extended-peak mode the slider value is the extended
        // peak; the operating point sits at the solved nominal Lw
        let lw = peakLuminance;
        if (extendedRange.extendedRangeActive && extendedRange.sliderSetsExtendedPeak) {
            const effectiveBlackLevel = blackLevelEnabled ? blackLevel : 0;
            lw = this.calculator.solveNominalFromExtendedPeak(peakLuminance, surroundLuminance, effectiveBlackLevel).lw;
        }

        const gamma = this.calculator.calculateSystemGamma(lw, surroundLuminance);
        this.current = { lw, ls: surroundLuminance, gamma };

        this.basicPath.attr("d", this.lineGenerator()(this.sampleCurveFrom(v => this.basicFormulaGamma(v, surroundLuminance))));
        this.extendedPath.attr("d", this.lineGenerator()(this.sampleCurveFrom(v => this.extendedRangeGamma(v, surroundLuminance))));
        this.implementedPath.attr("d", this.lineGenerator()(this.sampleImplementedCurve(surroundLuminance)));

        // Operating point (pinned to the axis edge if the solved Lw
        // falls outside it); nudge the label inside near the edges
        const cx = this.x(Math.min(this.LW_MAX, Math.max(this.LW_MIN, lw)));
        const xRange = this.x.range();
        this.operatingPoint.attr("transform", `translate(${cx},${this.y(gamma)})`);
        this.operatingLabel
            .attr("x", Math.min(Math.max(0, xRange[0] + 30 - cx), xRange[1] - 30 - cx))
            .call(t => drawLabel(t, `*γ* = ${gamma.toFixed(3)}`));

        // Nominal peak is this chart's x axis, so the dot's position
        // carries it; the surround is the fixed parameter every curve
        // here is drawn for, and is named nowhere else on the view.
        // The readout below reports the probe's peak, not the display's
        this.legendDetail.innerHTML =
            ` (nominal peak <var>L<sub>W</sub></var> = ${this.formatPeak(lw)} cd/m², ` +
            `${this.formatSurround(surroundLuminance)} cd/m² surround)`;

        // The probe persists across redraws, so it tracks the live
        // surround (setProbe re-clamps and keeps both controls in step)
        this.setProbe(this.probeLw);
    }

    /**
     * Format a peak luminance value for the hover readout
     */
    formatPeak(lw) {
        return Math.round(lw).toLocaleString("en-US").replace(/,/g, " ");
    }

    /**
     * Format a surround luminance value for the hover readout
     */
    formatSurround(ls) {
        if (ls >= 10) return ls.toFixed(0);
        if (ls >= 1) return ls.toFixed(1);
        return ls.toFixed(2);
    }

    /**
     * Render the crosshair, curve dots and readout at the probe
     * position. The probe persists and is re-rendered on every
     * redraw, so it follows the current surround luminance.
     */
    renderProbe() {
        if (!this.current) return;

        const lw = this.probeLw;
        const px = this.x(lw);
        const ls = this.current.ls;

        const implemented = this.calculator.calculateSystemGamma(lw, ls);
        const usesBasic = lw >= this.MONITORING_RANGE[0] && lw <= this.MONITORING_RANGE[1];
        const other = usesBasic ? this.extendedRangeGamma(lw, ls) : this.basicFormulaGamma(lw, ls);
        const otherName = usesBasic ? "Extended formula" : "Basic formula";

        this.hoverGroup.style("display", null);
        this.hoverLine.attr("x1", px).attr("x2", px);
        this.hoverDotImplemented.attr("transform", `translate(${px},${this.y(implemented)})`);
        this.hoverDotOther.attr("transform", `translate(${px},${this.y(other)})`);

        // "At", not "Nominal peak": that phrase names the display's own
        // setting, which the legend now carries, and this is the axis
        // position being read. Not "probe" either — in display work a
        // probe is a colorimeter, and nothing here is measured. The
        // surround comes from the legend too, so it is not repeated
        this.readout.innerHTML =
            `<span class="readout-item">At ${this.formatPeak(lw)} cd/m²</span> ` +
            `<span class="readout-item"><span class="readout-swatch swatch-adapted"></span>Implemented <var>γ</var> = ${implemented.toFixed(3)}</span> ` +
            `<span class="readout-item"><span class="readout-swatch swatch-reference"></span>${otherName} <var>γ</var> = ${other.toFixed(3)}</span>`;
    }
}


//=============================================================================
// GRAPH DISPLAY SETTINGS
//=============================================================================

/**
 * The calculator's display controls repeated in a graph view — the
 * Nominal Peak Luminance and Surround Luminance sliders with their
 * resets, the range mode's buttons and, in Extended Range, the pair that says which peak
 * the slider sets — so the settings a curve is drawn for can be changed
 * while the curve is in view.
 *
 * The repeated controls hold no state. A move writes the value to the
 * calculator's own slider and dispatches its "input" event, and a
 * button clicks the calculator's own, so the page updates by the single
 * path it always has. They re-read the originals on
 * hlg-values-updated — value, bounds, and the peak's label, which names
 * the extended peak in that mode — and the range buttons follow the
 * originals' clicks, which change the mode without emitting it.
 */
class GraphDisplaySettings {
    constructor(calculator, prefix) {
        this.calculator = calculator;

        // Every control of a view's block is the view's id prefix and
        // the control's role, so the block is addressed by prefix
        const own = role => document.getElementById(prefix + role);

        this.peak = {
            master: document.getElementById("peakLuminance"),
            masterLabel: document.getElementById("peakLuminanceTitle"),
            masterReset: document.getElementById("peakLuminanceResetBtn"),
            slider: own("DisplayPeak"),
            input: own("DisplayPeakInput"),
            label: own("DisplayPeakLabel"),
            reset: own("PeakResetBtn")
        };

        this.surround = {
            master: document.getElementById("surroundLuminance"),
            masterLabel: null,
            masterReset: document.getElementById("surroundLuminanceResetBtn"),
            slider: own("DisplaySurround"),
            input: own("DisplaySurroundInput"),
            label: null,
            reset: own("SurroundResetBtn")
        };

        this.pairs = [this.peak, this.surround];

        this.buttons = [
            { master: document.getElementById("nominalRangeButton"), button: own("NominalRangeBtn") },
            { master: document.getElementById("extendedRangeButton"), button: own("ExtendedRangeBtn") },
            { master: document.getElementById("peakNominalBtn"), button: own("PeakNominalBtn") },
            { master: document.getElementById("peakExtendedBtn"), button: own("PeakExtendedBtn") }
        ];

        // The peak pair belongs to Extended Range, so it is shown and
        // hidden with the calculator's own
        this.peakModeButtons = own("PeakModeButtons");
        this.peakModeToggle = document.getElementById("peakNominalOrExtendedToggle");
    }

    initialize() {
        this.pairs.forEach(pair => this.setupEventListeners(pair));

        // Registered after the ExtendedRangeController's, so the
        // original's own state is already set when the copy follows it
        this.buttons.forEach(({ master, button }) => {
            button.addEventListener("click", () => this.clickMaster(master, button));
            master.addEventListener("click", () => this.syncButtons());
        });

        document.addEventListener('hlg-values-updated', () => this.sync());
        this.sync();
    }

    setupEventListeners(pair) {
        pair.slider.addEventListener("input", () => this.setValue(pair, pair.slider.value));

        pair.reset.addEventListener("click", () => this.clickMaster(pair.masterReset, pair.reset));

        pair.input.addEventListener("input", () => {
            const typed = parseFloat(pair.input.value);
            if (Number.isFinite(typed)) this.setValue(pair, typed);
        });

        // A half-typed number was clamped on its way to the slider;
        // leaving the box shows what the slider took
        pair.input.addEventListener("blur", () => {
            pair.input.value = pair.master.value;
        });
    }

    /**
     * Set a value on the calculator's slider, which owns it, clamped
     * to its range as the calculator's own value box clamps
     */
    setValue(pair, value) {
        const min = parseFloat(pair.master.min);
        const max = parseFloat(pair.master.max);
        pair.master.value = Math.min(max, Math.max(min, parseFloat(value)));
        pair.master.dispatchEvent(new Event("input"));
    }

    /**
     * Re-read the calculator's sliders. A box being typed in keeps
     * what is typed; every other part follows.
     */
    sync() {
        this.pairs.forEach(pair => {
            pair.slider.min = pair.master.min;
            pair.slider.max = pair.master.max;
            pair.slider.step = pair.master.step;
            pair.slider.value = pair.master.value;
            pair.input.min = pair.master.min;
            pair.input.max = pair.master.max;
            if (document.activeElement !== pair.input) {
                pair.input.value = pair.master.value;
            }
            if (pair.label && pair.masterLabel) {
                pair.label.textContent = pair.masterLabel.textContent;
            }
        });
        this.syncButtons();
    }

    clickMaster(master, button) {
        this.holdPosition(button, () => master.click());
    }

    /**
     * Run a change that reaches the calculator, holding the view still.
     * The modes show and hide rows above this view — the 109% peak, its
     * toggle, the black level's slider and rows, the formulas and the
     * notes — so everything below them, this view included, would
     * otherwise move out from under the pointer by their height.
     * Scrolling by what the pressed button moved puts it back where it
     * was pressed.
     */
    holdPosition(button, change) {
        const before = button.getBoundingClientRect().top;
        change();
        window.scrollBy(0, button.getBoundingClientRect().top - before);
    }

    /**
     * Light the buttons the calculator's toggles have lit, and show the
     * peak pair when the calculator shows its own
     */
    syncButtons() {
        this.buttons.forEach(({ master, button }) => {
            button.classList.toggle("active", master.classList.contains("active"));
        });
        this.peakModeButtons.style.display =
            getComputedStyle(this.peakModeToggle).display === "none" ? "none" : "flex";
    }
}


//=============================================================================
// GRAPH PANEL CONTROLLER
//=============================================================================

/**
 * Owns the live graph panel: waits for d3 as HLGAdaptationCalculator
 * waits for KaTeX, initializes the views and switches between them.
 * makeGraphResponsive builds each chart once its container is
 * visible and has a width, hidden ones measuring 0; built views keep
 * redrawing while hidden, so switching needs no redraw. Without d3
 * the panel hides and the calculator works on.
 */
class GraphPanelController {
    constructor(calculator, uiController, extendedRangeController) {
        this.section = document.getElementById("graphPanelSection");

        // The signal chain over the views; its data-view lights the
        // blocks the current view is a graph of (see the CSS component)
        this.chainStrip = document.getElementById("chainStrip");

        this.cameraOetfGraph = new CameraOetfGraphController(calculator, uiController, extendedRangeController);
        this.eotfGraph = new EotfGraphController(calculator, uiController, extendedRangeController);
        this.surroundGammaGraph = new SurroundGammaGraphController(calculator, uiController, extendedRangeController);
        this.peakGammaGraph = new PeakGammaGraphController(calculator, uiController, extendedRangeController);

        // Each display-settings view's repeat of the display sliders;
        // they live with the panel, which hides as a whole without d3
        this.displaySettings = ["eotf", "surroundGamma", "peakGamma"]
            .map(prefix => new GraphDisplaySettings(calculator, prefix));

        // One entry per view, its toggle button and its container. The
        // EOTF opens the panel: it is the curve the calculator's own
        // settings shape (the others are the camera's and the two
        // system gamma relations)
        this.views = [
            { name: "cameraOetf", button: document.getElementById("graphViewCameraOetfBtn"), element: document.getElementById("cameraOetfView") },
            { name: "eotf", button: document.getElementById("graphViewEotfBtn"), element: document.getElementById("eotfView") },
            { name: "surroundGamma", button: document.getElementById("graphViewSurroundGammaBtn"), element: document.getElementById("surroundGammaView") },
            { name: "peakGamma", button: document.getElementById("graphViewPeakGammaBtn"), element: document.getElementById("peakGammaView") }
        ];
    }

    /**
     * Initialize the panel once d3 is available
     */
    initialize(attempts = 0) {
        const maxAttempts = 20; // 20 * 50ms = 1 second max wait

        if (window.d3) {
            this.cameraOetfGraph.initialize();
            this.eotfGraph.initialize();
            this.surroundGammaGraph.initialize();
            this.peakGammaGraph.initialize();
            this.displaySettings.forEach(settings => settings.initialize());
            this.setupViewToggle();
        } else if (attempts < maxAttempts) {
            setTimeout(() => this.initialize(attempts + 1), 50);
        } else {
            console.warn("d3 failed to load after multiple attempts. The graph panel is hidden.");
            this.section.style.display = "none";
        }
    }

    /**
     * Set up the view toggle buttons
     */
    setupViewToggle() {
        this.views.forEach(view => {
            view.button.addEventListener("click", () => this.setView(view.name));
        });
    }

    /**
     * Show one view and hide the others, and light the blocks of
     * the signal chain the shown view is a graph of
     */
    setView(name) {
        this.views.forEach(view => {
            view.button.classList.toggle("active", view.name === name);
            view.element.style.display = view.name === name ? "block" : "none";
        });
        this.chainStrip.dataset.view = name;
    }
}
