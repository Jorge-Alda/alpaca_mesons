"""Extract literal tables from the cached arXiv source archives."""
from pathlib import Path
import re, tarfile, json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/odd_parity"
OUT.mkdir(parents=True, exist_ok=True)

def source(key, filename):
    with tarfile.open(ROOT / "memory/sources" / key / "arxiv-source.tar") as tf:
        return tf.extractfile(filename).read().decode("latin1")

def cached_text(relative):
    return (ROOT / relative).read_text(encoding="utf-8")

def require_once(pattern, text, description, flags=0):
    matches = list(re.finditer(pattern, text, flags))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {description}; found {len(matches)}")
    return matches[0]

def parse_kloe_fit(text):
    label = require_once(r"\\label\{param\}", text, "KLOE fit-results label")
    start = text.rfind("\\begin{table}", 0, label.start())
    end = text.find("\\end{table}", label.end())
    if start < 0 or end < 0:
        raise RuntimeError("KLOE fit-results table boundaries not found")
    table = text[start:end]
    if not re.search(r"\\caption\{Fit results\.\s*Masses and widths are in MeV\.", table):
        raise RuntimeError("unexpected KLOE fit-results caption")
    if not re.search(r"parameter\s*&\s*fit\(a\)\s*&\s*fit\(b\)\s*&\s*fit\(c\)", table):
        raise RuntimeError("unexpected KLOE fit-results header")
    names = {
        r"M_{\rho^0}": ("rho0_mass", "MeV"),
        r"M_{\rho^+}": ("rho_plus_mass", "MeV"),
        r"M_{\rho^-}": ("rho_minus_mass", "MeV"),
        r"\Gamma_{\rho^0}": ("rho0_width", "MeV"),
        r"\Gamma_{\rho^+}": ("rho_plus_width", "MeV"),
        r"\Gamma_{\rho^-}": ("rho_minus_width", "MeV"),
        r"a_d": ("direct_amplitude", "dimensionless"),
        r"\phi_d": ("direct_phase", "rad"),
        r"a_{\omega}\times10^3": ("omega_pi_amplitude_x1e3", "dimensionless"),
        r"\phi_{\omega}": ("omega_pi_phase", "rad"),
    }
    expected_order = [item[0] for item in names.values()]
    expected_blank = {"rho0_mass", "rho_minus_mass", "rho0_width", "rho_minus_width"}
    parsed = {}
    tabular_start = table.find("\\begin{tabular}")
    tabular_end = table.find("\\end{tabular}", tabular_start)
    if tabular_start < 0 or tabular_end < 0:
        raise RuntimeError("KLOE fit-results tabular block not found")
    for line in table[tabular_start:tabular_end].splitlines():
        if "&" not in line or "\\\\" not in line:
            continue
        fields = [re.sub(r"\\(?:hline|cline\{[^}]*\})", "", item).strip() for item in line.split("&")]
        parameter = fields[0].replace("$", "")
        if parameter in {"parameter", r"\chi^2 [p(\chi^2)]"}:
            continue
        if parameter not in names:
            raise RuntimeError(f"unexpected KLOE fit parameter row: {parameter!r}")
        if len(fields) != 4:
            raise RuntimeError(f"malformed KLOE fit row: {line!r}")
        key, unit = names[parameter]
        if key in parsed:
            raise RuntimeError(f"duplicate KLOE fit parameter row: {key}")
        cell = fields[1].replace("$", "").replace("\\\\", "").strip()
        if not cell:
            if key not in expected_blank:
                raise RuntimeError(f"unexpected blank KLOE fit(a) cell: {key}")
            parsed[key] = [key, "", "", "", unit]
            continue
        match = re.fullmatch(r"([+-]?\d+(?:\.\d*)?)\s*\\pm\s*(\d+(?:\.\d*)?)\s*\\pm\s*(\d+(?:\.\d*)?)", cell)
        if not match or key in expected_blank:
            raise RuntimeError(f"unexpected KLOE fit(a) value: {key}={cell!r}")
        parsed[key] = [key, match.group(1), match.group(2), match.group(3), unit]
    if list(parsed) != expected_order:
        raise RuntimeError(f"KLOE fit parameter rows changed: {list(parsed)!r}")
    return list(parsed.values())

def parse_pdg_omega_branching(text):
    heading = require_once(r"ω\(782\)\s+DECAY MODES", text, "PDG omega decay-mode heading")
    end = text.find("Γ3", heading.end())
    if end < 0:
        raise RuntimeError("PDG omega decay-mode boundary Γ3 not found")
    section = text[heading.end():end]
    match = require_once(r"Γ2\s+π0\s+γ\s+\(\s*([0-9.]+)\s*±\s*([0-9.]+)\s*\)\s*%", section,
                         "PDG omega-to-pi0-gamma branching datum")
    return f"{float(match.group(1)) / 100:g}", f"{float(match.group(2)) / 100:g}"

def write(name, text, meta):
    (OUT / name).write_text(text, encoding="utf-8")
    (OUT / (Path(name).stem + ".metadata.json")).write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

def err(s):
    a = re.findall(r"([0-9.]+)\\s*\\pm", s)
    return a

def likelihood(diagonal_sources, nuisance_sources=()):
    """Versioned, user-authorized first-release likelihood metadata."""
    excluded = [item["source_id"] for item in nuisance_sources]
    return {
        "schema_version": "odd-parity-likelihood-v1",
        "family": "diagonal_gaussian",
        "authorization": "User-authorized diagonal published errors where covariance is absent.",
        "asymmetric_error_rule": "mean_absolute",
        "diagonal_source_ids": list(diagonal_sources),
        "nuisance_sources": list(nuisance_sources),
        "nuisance_source_ids_excluded_from_diagonal": excluded,
    }

# CLEO: retain the published Qbar^2 F points, with separate statistical and systematic errors.
cleo = source("CLEO:1997fho", "clns_97-1477.tex")
rows = []
for label, meson in [("tab:piz_table_1", "pi0"),
                     ("tab:etaz_table_1", "eta"), ("tab:etaz_table_2", "eta"), ("tab:etaz_table_3", "eta"),
                     ("tab:etap_table_1", "eta_prime"), ("tab:etap_table_2", "eta_prime"),
                     ("tab:etap_table_3", "eta_prime"), ("tab:etap_table_4", "eta_prime"),
                     ("tab:etap_table_5", "eta_prime"), ("tab:etap_table_6", "eta_prime")]:
    block = cleo[cleo.index("\\label{" + label + "}"):]
    block = block[:block.index("\\end{table}")]
    for line in block.splitlines():
        if "&" not in line or "\\\\" not in line or "\\pm" not in line:
            continue
        fields = [x.replace("$", "").replace("{", "").replace("}", "").strip() for x in line.split("&")]
        if len(fields) < 8:
            continue
        q = re.search(r"([0-9.]+)", fields[5])
        f = re.findall(r"[0-9.]+", fields[7])
        if q and len(f) >= 3:
            rows.append([meson, label, q.group(1),
                         f"{float(f[0]) * 0.01:g}", f"{float(f[1]) * 0.01:g}", f"{float(f[2]) * 0.01:g}"])
write("cleo_tff.csv", "meson,source_table,q2_GeV2,q2F_GeV,stat_GeV,syst_GeV\n" + "\n".join(",".join(r) for r in rows) + "\n", {
    "bibtex_key": "CLEO:1997fho", "source_location": "Source labels tab:piz_table_1, tab:etaz_table_1--3, and tab:etap_table_1--6 in clns_97-1477.tex, pp. 24-26",
    "observable": "Qbar^2 |F(Qbar^2)|", "units": {"q2": "GeV^2", "q2F": "GeV"},
    "source_scale": "The source table prints q2F and its errors in 0.01 GeV; CSV values are converted to GeV by multiplying by 0.01.",
    "errors": "stat_GeV and syst_GeV are the separately quoted diagonal errors after the same 0.01 GeV conversion; no covariance matrix supplied; normalization/luminosity systematics are not reconstructed as bin correlations.",
    "likelihood": likelihood(["cleo-statistical", "cleo-point-systematic"])
})

# BaBar Table 11, including the quoted combined TFF error and the two cross-section components.
babar = source("BaBar:2009rrj", "babar-pub-09006.tex")
block = babar[babar.index("\\label{tab11}"):]
block = block[:block.index("\\end{table}")]
rows = []
for line in block.splitlines():
    if "&" not in line or "\\\\" not in line or "\\pm" not in line:
        continue
    f = [x.replace("$", "").replace("{", "").replace("}", "").strip() for x in line.split("&")]
    if len(f) != 4:
        continue
    nums = re.findall(r"([0-9.]+)\s*\\pm", f[2])
    tff = re.search(r"([0-9.]+)[.]?\s*\\pm\s*([0-9.]+)", f[3])
    asym = re.search(r"([0-9.]+)[.]?_(-?[0-9.]+)\.\^\+([0-9.]+)", f[3])
    if len(nums) == 2 and (tff or asym):
        if tff: rows.append([f[0], f[1], nums[0], nums[1], tff.group(1), tff.group(2), tff.group(2)])
        else: rows.append([f[0], f[1], nums[0], nums[1], asym.group(1), str(abs(float(asym.group(2)))), asym.group(3)])
write("babar_pi0_table11.csv", "q2_interval_GeV2,q2bar_GeV2,dsigma_stat_fb_per_GeV2,dsigma_syst_fb_per_GeV2,q2F_MeV,q2F_err_low_MeV,q2F_err_high_MeV\n" + "\n".join(",".join(r) for r in rows) + "\n", {
    "bibtex_key": "BaBar:2009rrj", "source_location": "Table 11, arXiv source babar-pub-09006.tex, p. 19",
    "units": {"q2": "GeV^2", "cross_section": "fb/GeV^2", "q2F": "MeV"},
    "errors": "Cross section has statistical and Q2-dependent systematic columns; q2F has combined error. For an asymmetric q2F entry, the eventual authorized diagonal sigma is (q2F_err_low_MeV + q2F_err_high_MeV)/2 while the raw lower/upper magnitudes remain separate. The paper quotes additional Q2-independent 3% cross-section and 2.3% form-factor normalizations.",
    "likelihood": likelihood(
        ["babar-2009-point-error"],
        [
            {"source_id": "babar-2009-cross-section-normalization", "fraction": 0.03,
             "response": "multiplicative", "predictor_connection": "absent"},
            {"source_id": "babar-2009-form-factor-normalization", "fraction": 0.023,
             "response": "multiplicative", "predictor_connection": "absent"},
        ],
    )
})

# BaBar ISR spectrum: table is two side-by-side blocks; parse each literal row.
isr = source("BaBar:2004ytv", "paper.tex")
block = isr[isr.index("\\label{sumtab}"):]
block = block[:block.index("\\end{table*}")]
rows = []
def isr_side(fields):
    if len(fields) != 5: return None
    mass = re.search(r"([0-9.]+)", fields[0])
    n = re.search(r"(-?[0-9]+)\s*\\pm\s*([0-9]+)\s*\\pm\s*([0-9]+)", fields[1])
    eff = re.search(r"([0-9.]+)\s*\\pm\s*([0-9.]+)", fields[2])
    lum = re.search(r"([0-9]+)", fields[3])
    sig = re.search(r"(-?[0-9.]+)\s*\\pm\s*([0-9.]+)\s*\\pm\s*([0-9.]+)", fields[4])
    if not all((mass, n, eff, lum, sig)): return None
    return [mass.group(1), n.group(1), n.group(2), n.group(3), eff.group(1), eff.group(2), lum.group(1),
            sig.group(1), sig.group(2), sig.group(3)]
for line in block.splitlines():
    if "&" not in line: continue
    fields = [x.replace("$", "").strip() for x in line.split("&")]
    if len(fields) != 10: continue
    left, right = isr_side(fields[:5]), isr_side(fields[5:])
    if left: rows.append(left)
    if right: rows.append(right)
write("babar_3pi_born_spectrum.csv", "sqrt_s_GeV,Ncorr,Ncorr_stat,Ncorr_syst,efficiency_percent,efficiency_err_percent,dL_nb_inv,sigma_nb,sigma_stat_nb,sigma_syst_nb\n" + "\n".join(",".join(r) for r in rows) + "\n", {
    "bibtex_key": "BaBar:2004ytv", "source_location": "Table sumtab, arXiv source paper.tex, pp. 25-26",
    "units": {"sqrt_s": "GeV", "efficiency": "%", "efficiency_error": "%", "dL": "nb^-1", "sigma": "nb"},
    "errors": "The CSV retains Ncorr statistical/systematic errors, efficiency and its quoted uncertainty, and sigma statistical/systematic errors. For sigma, the source assigns the systematic component to event counts, detection efficiencies, and radiative correction, and states that mass-bin systematic errors are fully correlated. The 1.2% luminosity nuisance is retained separately; it is quoted explicitly in the fitted-product uncertainty discussion, not folded into the per-bin sigma systematic column.",
    "covariance_assumption": "Only the published sigma statistical error is diagonal. The published sigma systematic vector is fully correlated across mass bins and is represented by one additive standardized nuisance; luminosity is a separate multiplicative nuisance.",
    "likelihood": likelihood(
        ["babar-2004-statistical"],
        [
            {"source_id": "babar-2004-mass-bin-systematic", "response": "additive",
             "response_column": "sigma_syst_nb", "unit": "nb", "predictor_connection": "absent"},
            {"source_id": "babar-2004-luminosity-normalization", "fraction": 0.012,
             "response": "multiplicative", "predictor_connection": "absent"},
        ],
    )
})

# KLOE publishes fit parameters, not the 1874 bin contents; preserve the tabulated fit-a inputs.
kloe = source("KLOE:2003kas", "articolo_errata.tex")
kloe_rows = parse_kloe_fit(kloe)
write("kloe_phi_3pi_dalitz_fit.csv", "parameter,value,stat_error,syst_error,unit\n" + "\n".join(",".join(row) for row in kloe_rows) + "\n", {
    "bibtex_key": "KLOE:2003kas", "source_location": "Table (fit results), arXiv source articolo_errata.tex, pp. 7-8",
    "units": "Per unit column", "errors": "Fit statistical and systematic errors as printed. Blank rho0/rho- mass and width cells are preserved exactly from fit(a), whose isospin/CPT constraints leave the fitted value printed only on the rho+ row; they are not missing measurements.",
    "likelihood": likelihood(["kloe-fit-statistical", "kloe-fit-systematic"]),
    "fit_readiness": "No covariance matrix or bin-level Dalitz data supplied; this is a fit-parameter table, not a reconstructed Dalitz sample. CMD-2 tabulated Dalitz inputs are unavailable in the cached Tier-0 material."
})

pdg = cached_text("memory/sources/ParticleDataGroup:2024cfk/rpp2024-list-omega-782.txt")
pdg_value, pdg_error = parse_pdg_omega_branching(pdg)
write("omega_pi0_gamma_pdg.csv", "observable,value,stat_or_total_error,unit\n" +
      f"omega_to_pi0_gamma_branching_fraction,{pdg_value},{pdg_error},fraction\n", {
    "bibtex_key": "ParticleDataGroup:2024cfk", "source_location": "PDG 2024 omega(782) listing, p. 3, decay modes",
    "units": {"branching_fraction": "fraction"}, "source_scale": "PDG prints (8.35 +/- 0.27) percent; CSV converts the percentage and error to fraction units by dividing by 100.", "errors": "PDG quoted uncertainty; this is an averaged branching fraction.",
    "likelihood": likelihood(["pdg-2024-omega-pi0-gamma-total"])
})
