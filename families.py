"""Outcome-family regexes, placebo-group rules, gating and analysis-method rules (fixed source-project family definitions)."""
import re

I = re.IGNORECASE

# family -> (category, improvement_sign, title_regex, eligibility_regex)
# improvement_sign: +1 means a decrease is improvement (response = -change); -1 means increase is improvement.
FAMILIES = {
    "MADRS":   ("CR", +1, r"\bMADRS\b|Montgomery[\s-]*(?:and[\s-]*)?[AÅ]sberg", r"\bMADRS\b|Montgomery[\s-]*(?:and[\s-]*)?[AÅ]sberg"),
    "HAMD":    ("CR", +1, r"\bHAM[\s-]?D(?:-?\d+)?\b|\bHDRS\b|\bHRSD\b|Hamilton\s+(?:Rating\s+Scale\s+for\s+|Depression)", r"\bHAM[\s-]?D(?:-?\d+)?\b|\bHDRS\b|\bHRSD\b|Hamilton\s+(?:Rating\s+Scale\s+for\s+)?Depression"),
    "PANSS":   ("CR", +1, r"\bPANSS\b(?=.*\btotal\b)|(?=.*\btotal\b).*\bPANSS\b|Positive\s+and\s+Negative\s+Syndrome\s+Scale(?=.*total)", r"\bPANSS\b|Positive\s+and\s+Negative\s+Syndrome"),
    "HAMA":    ("CR", +1, r"\bHAM[\s-]?A\b|Hamilton\s+(?:Rating\s+Scale\s+for\s+)?Anxiety", r"\bHAM[\s-]?A\b|Hamilton\s+(?:Rating\s+Scale\s+for\s+)?Anxiety"),
    "YMRS":    ("CR", +1, r"\bYMRS\b|Young\s+Mania", r"\bYMRS\b|Young\s+Mania"),
    "ADHDRS":  ("CR", +1, r"\bADHD[\s-]*RS\b|ADHD\s+Rating\s+Scale", r"\bADHD[\s-]*RS\b|ADHD\s+Rating\s+Scale"),
    "PAIN":    ("PR", +1, r"(?=.*\bpain\b)(?=.*(?:\bNRS\b|numeric(?:al)?\s+rating|\bVAS\b|visual\s+analog|pain\s+intensity|pain\s+score|average\s+(?:daily\s+)?pain))(?!.*WOMAC)(?!.*(?:BPI\s+interference|interference))", r"(?:\bpain\b).{0,60}(?:\bNRS\b|numeric|\bVAS\b|visual\s+analog|score|intensity|rating)|(?:\bNRS\b|numeric\s+rating|\bVAS\b).{0,60}\bpain\b"),
    "WOMACP":  ("PR", +1, r"WOMAC(?=.*\bpain\b)", r"WOMAC"),
    "IPSS":    ("PR", +1, r"\bIPSS\b|International\s+Prostate\s+Symptom", r"\bIPSS\b|International\s+Prostate\s+Symptom"),
    "HOTFLASH":("PR", +1, r"hot\s+fl(?:a|u)sh|vasomotor\s+symptom", r"hot\s+fl(?:a|u)sh|vasomotor"),
    "HBA1C":   ("OB", +1, r"\bHbA1c\b|\bA1C\b|glyc(?:at|osyl)ated\s+ha?emoglobin|ha?emoglobin\s+A1c", r"\bHbA1c\b|\bA1C\b|glyc(?:at|osyl)ated\s+ha?emoglobin|ha?emoglobin\s+A1c"),
    "FPG":     ("OB", +1, r"\bFPG\b|fasting\s+(?:plasma\s+|serum\s+|blood\s+)?glucose", r"\bFPG\b|fasting\s+(?:plasma\s+|serum\s+|blood\s+)?glucose"),
    "LDL":     ("OB", +1, r"\bLDL\b(?:[\s-]*C\b|[\s-]*cholesterol)|low[\s-]density\s+lipoprotein", r"\bLDL\b|low[\s-]density\s+lipoprotein"),
    "FEV1":    ("OB", -1, r"\bFEV\s*1\b|\bFEV1\b|forced\s+expiratory\s+volume\s+in\s+(?:1|one)", r"\bFEV\s*1\b|\bFEV1\b|forced\s+expiratory\s+volume"),
    "SBP":     ("OBB", +1, r"(?=.*(?:systolic\s+blood\s+pressure|\bSBP\b|\bSeSBP\b|\bmsSBP\b|\bSiSBP\b))(?!.*(?:ambulatory|\bABPM\b|24[\s-]*h(?:ou)?r|daytime|night[\s-]*time))", r"systolic\s+blood\s+pressure|\bSBP\b|\bmsSBP\b|\bSeSBP\b|\bSiSBP\b"),
    "SBPABPM": ("OBB", +1, r"(?=.*(?:systolic|\bSBP\b))(?=.*(?:ambulatory|\bABPM\b|24[\s-]*h(?:ou)?r))", r"systolic|\bSBP\b"),
    "WEIGHT":  ("OBB", +1, r"\bbody\s+weight\b(?!.*(?:percent|%|kg/m))|\bweight\s+\(kg\)", r"\bbody\s+weight\b|\bBMI\b|body\s+mass\s+index"),
}

FAM_TITLE = {k: re.compile(v[2], I) for k, v in FAMILIES.items()}
FAM_ELIG = {k: re.compile(v[3], I) for k, v in FAMILIES.items()}

# Families in priority order when a title matches several (more specific first).
FAMILY_ORDER = ["SBPABPM", "SBP", "WOMACP", "PANSS", "MADRS", "HAMA", "HAMD", "YMRS", "ADHDRS", "IPSS",
                "HOTFLASH", "PAIN", "HBA1C", "FPG", "LDL", "FEV1", "WEIGHT"]

CHANGE_WORD = re.compile(r"change|reduction|decrease|improvement|difference\s+from\s+baseline", I)
PERCENT_WORD = re.compile(r"percent(?:age)?\s+change|%\s*change|percent\s+reduction", I)
PARAM_OK = {"MEAN", "LEAST_SQUARES_MEAN"}

PLACEBO = re.compile(r"placebo|vehicle|\bsham\b|dummy", I)
COMBO = re.compile(r"\+|\bplus\b|add[\s-]?on|adjunct|background|in\s+combination|\bwith\s", I)

THRESH = re.compile(r"(?:≥|>=|=>|>|≤|<=|=<|<|at\s+least|greater\s+than|more\s+than|less\s+than|higher\s+than|lower\s+than|between|minimum\s+of|\bof\s+\d|\d+(?:\.\d+)?\s*(?:to|-|–)\s*\d+)", I)
EXCL_SPLIT = re.compile(r"exclusion\s+criteria", I)

M_LOCF = re.compile(r"\bLOCF\b|last\s+observation\s+carried", I)
M_MMRM = re.compile(r"\bMMRM\b|mixed[\s-]*(?:effects?\s+)?model|repeated[\s-]+measures|linear\s+mixed", I)
M_MI = re.compile(r"multiple\s+imputation|\bMI\b", I)
M_OC = re.compile(r"observed\s+cases?|\bOC\b|completers?", I)


def classify_title(title):
    for fam in FAMILY_ORDER:
        if FAM_TITLE[fam].search(title or ""):
            return fam
    return None


def gated_v1(elig_text, fam, window=120):
    """v1 rule (superseded by v1.1 `gated`)."""
    if not elig_text:
        return 0
    inc = EXCL_SPLIT.split(elig_text, maxsplit=1)[0]
    for m in FAM_ELIG[fam].finditer(inc):
        lo, hi = max(0, m.start() - window), min(len(inc), m.end() + window)
        seg = inc[lo:hi]
        if THRESH.search(seg) and re.search(r"\d", seg):
            return 1
    return 0


# ---- v1.1 amendments (PREREG_placebo_v1_1_amendment.md) ----
# A1: gating on the same CONSTRUCT (any scale of the same condition domain), not only the same scale.
CONSTRUCT = {
    "MADRS": "DEP", "HAMD": "DEP", "HAMA": "ANX", "PANSS": "PSY", "YMRS": "MANIA", "ADHDRS": "ADHD",
    "PAIN": "PAIN", "WOMACP": "PAIN", "IPSS": "IPSS", "HOTFLASH": "VMS", "HBA1C": "GLY", "FPG": "GLY",
    "LDL": "LIPID", "FEV1": "LUNG", "SBP": "BP", "SBPABPM": "BP", "WEIGHT": "WT",
}
CONSTRUCT_ELIG = {
    "DEP": r"\bMADRS\b|Montgomery|\bHAM[\s-]?D|\bHDRS\b|\bHRSD\b|Hamilton\s+(?:Rating\s+Scale\s+for\s+)?Depression|\bIDS(?:-SR|-C)?\b|\bQIDS\b|\bCGI[\s-]*S\b|Clinical\s+Global\s+Impression",
    "ANX": r"\bHAM[\s-]?A\b|Hamilton\s+(?:Rating\s+Scale\s+for\s+)?Anxiety|\bCGI[\s-]*S\b|Clinical\s+Global\s+Impression",
    "PSY": r"\bPANSS\b|Positive\s+and\s+Negative\s+Syndrome|\bBPRS\b|\bCGI[\s-]*S\b|Clinical\s+Global\s+Impression",
    "MANIA": r"\bYMRS\b|Young\s+Mania|\bCGI[\s-]*(?:BP|S)\b|Clinical\s+Global\s+Impression",
    "ADHD": r"\bADHD[\s-]*RS\b|ADHD\s+Rating\s+Scale|\bCAARS\b|\bCGI[\s-]*S\b|Clinical\s+Global\s+Impression",
    "PAIN": r"\bpain\b|\bNRS\b|numeric(?:al)?\s+rating|\bVAS\b|visual\s+analog|WOMAC|\bBPI\b",
    "IPSS": r"\bIPSS\b|International\s+Prostate\s+Symptom|\bAUA[\s-]*SI\b",
    "VMS": r"hot\s+fl(?:a|u)sh|vasomotor",
    "GLY": r"\bHbA1c\b|\bA1C\b|glyc(?:at|osyl)ated\s+ha?emoglobin|ha?emoglobin\s+A1c|\bFPG\b|fasting\s+(?:plasma\s+|serum\s+|blood\s+)?glucose",
    "LIPID": r"\bLDL|low[\s-]density\s+lipoprotein|\bcholesterol\b|\btriglycerides?\b",
    "LUNG": r"\bFEV\s*1|\bFEV1\b|forced\s+expiratory\s+volume|\bPEF\b|peak\s+expiratory",
    "BP": r"systolic|diastolic|blood\s+pressure|\bSBP\b|\bDBP\b|\bmsSBP\b|\bSeSBP\b|\bSiSBP\b|\bmsDBP\b|\bSeDBP\b",
    "WT": r"\bbody\s+weight\b|\bBMI\b|body\s+mass\s+index",
}
CONSTRUCT_RE = {k: re.compile(v, I) for k, v in CONSTRUCT_ELIG.items()}
# A4: numeric threshold = comparator immediately followed by a number, or a numeric range,
#     and the number is not a time duration / age.
NUM_THRESH = re.compile(
    r"(?:≥|>=|=>|\\?>|≤|<=|=<|\\?<|at\s+least|(?:greater|more|less|higher|lower)\s+than(?:\s+or\s+equal\s+to)?|(?:greater|more|less|higher|lower)\s+(?:than\s+)?or\s+equal\s+to|"
    r"between|minimum(?:\s+of)?|maximum(?:\s+of)?|of|score\s+of|exceed(?:ing|s)?|above|below)\s*(?:\(?\\?[<>]=?\)?\s*)?(\d+(?:\.\d+)?)"
    r"|(\d+(?:\.\d+)?)\s*(?:%|mg/dl|mmol/l|mmhg|mm\s*hg|kg/m)?\s*(?:to|-|–|and)\s*\d+(?:\.\d+)?", I)
NOT_GATE_UNIT = re.compile(r"^\s*(?:\)|\]|%\s*of\s+days)?\s*(?:years?|yrs?|months?|weeks?|days?|hours?|hrs?|times?|cigarettes?|pack|episodes?\s+per\s+(?:year|month)|visits?|doses?|ppm)\b", I)
BULLET_SPLIT = re.compile(r"(?:\n|\s\*\s|•|\s\d{1,2}[.)]\s)")


def gated(elig_text, fam):
    """v1.1: inclusion criteria contain, within the same bullet, a construct term and a numeric threshold."""
    if not elig_text:
        return 0
    inc = EXCL_SPLIT.split(elig_text, maxsplit=1)[0]
    cre = CONSTRUCT_RE[CONSTRUCT[fam]]
    for bullet in BULLET_SPLIT.split(inc):
        if not cre.search(bullet):
            continue
        for m in NUM_THRESH.finditer(bullet):
            tail = bullet[m.end():m.end() + 25]
            if NOT_GATE_UNIT.search(tail):
                continue
            # age criteria
            head = bullet[max(0, m.start() - 25):m.start()]
            if re.search(r"\bage[ds]?\b|years?\s+old", head + tail, I):
                continue
            return 1
    return 0


# A2: weight-gain-is-improvement context (cachexia, anorexia, wasting, malnutrition, underweight)
WEIGHT_GAIN_CTX = re.compile(r"cachexia|anorexia|wasting|malnutrition|underweight|sarcopenia|weight\s+loss\s+(?:due|associated|in\s+patients)|involuntary\s+weight\s+loss|failure\s+to\s+thrive", I)


def analysis_method(texts):
    t = " ".join(x for x in texts if x)
    if M_LOCF.search(t):
        return "LOCF"
    if M_MMRM.search(t):
        return "MMRM"
    if M_MI.search(t):
        return "MI"
    if M_OC.search(t):
        return "OC"
    return "UNKNOWN"
