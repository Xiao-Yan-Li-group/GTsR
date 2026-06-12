import base64
import hashlib
import logging
import sys
import tempfile
from io import BytesIO
from pathlib import Path

import streamlit as st
from ase.build import make_supercell
from ase.io import read, write
from stmol import py3Dmol, showmol

from streamlit_extras.colored_header import colored_header

WEBAPP_DIR = Path(__file__).resolve().parents[2]
PROJECT_DIR = WEBAPP_DIR.parent
MAX_UPLOAD_BYTES = 10 * 1024 * 1024
VIEWER_WIDTH = 800
VIEWER_HEIGHT = 500
LOGGER = logging.getLogger(__name__)

logs_dir = WEBAPP_DIR / "logs"
logs_dir.mkdir(parents=True, exist_ok=True)

if "gtsr_workspace" not in st.session_state:
    st.session_state["gtsr_workspace"] = tempfile.TemporaryDirectory(
        prefix="session-",
        dir=logs_dir,
    )

session_dir = Path(st.session_state["gtsr_workspace"].name)
input_dir = session_dir / "uploads"
output_dir = session_dir / "predictions"
input_dir.mkdir(parents=True, exist_ok=True)
output_dir.mkdir(parents=True, exist_ok=True)

if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from gtsr import GTsRunner

@st.cache_resource(show_spinner=False)
def get_runner(checkpoint: str) -> GTsRunner:
    return GTsRunner(checkpoint=checkpoint)


def show_high_res_image(image_path: Path, display_width: int = 600) -> None:
    encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
    st.html(
        f'<img src="data:image/png;base64,{encoded}" alt="GTsR logo" '
        f'style="width:{display_width}px;max-width:100%;height:auto;">'
    )


@st.cache_data(show_spinner=False)
def read_cif(cif_data: bytes):
    return read(BytesIO(cif_data), format="cif")


@st.cache_data(show_spinner=False)
def normalize_cif(cif_data: bytes) -> str:
    structure = read_cif(cif_data)
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    normalized = make_supercell(structure, identity)

    cif_buffer = BytesIO()
    write(cif_buffer, normalized, format="cif")
    return cif_buffer.getvalue().decode("utf-8")


def cif_viewer(
    cif_text: str,
    style: str,
    background: str,
    show_labels: bool,
    show_cell: bool,
    spin: bool,
):
    viewer = py3Dmol.view(width=VIEWER_WIDTH, height=VIEWER_HEIGHT)
    viewer.addModel(cif_text, "cif")

    if style == "Ball and stick":
        viewer.setStyle({"stick": {"radius": 0.16}, "sphere": {"scale": 0.28}})
    elif style == "Space filling":
        viewer.setStyle({"sphere": {"scale": 0.75}})
    elif style == "Line":
        viewer.setStyle({"line": {"linewidth": 1.5}})
    else:
        viewer.setStyle({"stick": {"radius": 0.18}})

    if show_labels:
        viewer.addPropertyLabels(
            "elem",
            {},
            {
                "fontColor": "#003D7C",
                "fontSize": 12,
                "backgroundColor": "#FFF9E8",
                "backgroundOpacity": 0.75,
            },
        )

    viewer.setBackgroundColor(background)
    if show_cell:
        viewer.addUnitCell()

    viewer.zoomTo()

    if spin:
        viewer.spin("y", 0.6)

    return viewer


def show_structure_metrics(structure, key: str) -> None:
    with st.container(key=key):
        formula_col, atoms_col = st.columns(2)
        formula_col.metric("Chemical Formula", structure.get_chemical_formula())
        atoms_col.metric("Number of Atoms", len(structure))

    st.html(
        f"""
        <style>
        .st-key-{key} [data-testid="stMetricValue"] {{
            font-size: 1.25rem;
        }}
        </style>
        """
    )


def show_cleaned_structure(
    cif_path: str,
    style: str,
    show_labels: bool,
    spin: bool,
    metric_key: str,
) -> None:
    path = Path(cif_path)
    cleaned_cif = path.read_text(encoding="utf-8")
    structure = read(path)
    show_structure_metrics(structure, metric_key)
    cleaned_viewer = cif_viewer(
        cif_text=cleaned_cif,
        style=style,
        background="#FFFFFF",
        show_labels=show_labels,
        show_cell=True,
        spin=spin,
    )
    showmol(cleaned_viewer, height=VIEWER_HEIGHT, width=VIEWER_WIDTH)


st.title("GTSR")
st.html(
    '<p style="font-size:1rem;line-height:1.8;margin:0 0 1rem 0;">'
    'A <strong style="color:#D71920;font-size:1.5em;line-height:0;vertical-align:baseline;">G</strong>NN Based '
    '<strong style="color:#D71920;font-size:1.5em;line-height:0;vertical-align:baseline;">T</strong>ool for '
    '<strong style="color:#D71920;font-size:1.5em;line-height:0;vertical-align:baseline;">S</strong>olvent '
    '<strong style="color:#D71920;font-size:1.5em;line-height:0;vertical-align:baseline;">R</strong>emoval from MOF with stability check.'
    "</p>"
)

show_high_res_image(WEBAPP_DIR / "imgs" / "gtsr_logo.png", display_width=600)

uploaded_file = st.file_uploader(label="Upload Your **C**rystallographic **I**nformation **F**ile", type=["cif"], help="Upload Raw Structure",)

if uploaded_file is not None:

    cif_data = uploaded_file.getvalue()

    with st.spinner("Reading CIF..."):
        structure = read_cif(cif_data)
        cif_text = normalize_cif(cif_data)

    colored_header(
                label="Raw Structure",
                description="",
                color_name="orange-70",
            )

    type_col, label_col, spin_col = st.columns(3)

    with type_col:
        render_style = st.selectbox(
                                    "Style",
                                    ["Ball and stick", "Stick", "Space filling", "Line"],
                                    key="raw_rs"
                                )
    with label_col:
        show_labels = st.selectbox(
                                    "Element",
                                    [False, True],
                                    format_func=lambda value: "Show" if value else "Hide",
                                    key="raw_sl"
                                )
    with spin_col:
        spin = st.selectbox(
                            "Rotation",
                            [False, True],
                            format_func=lambda value: "On" if value else "Off",
                            key="raw_sp"
                                )

    show_structure_metrics(structure, "strucinfo_raw")

    viewer = cif_viewer(
        cif_text=cif_text,
        style=render_style,
        background="#FFFFFF",
        show_labels=show_labels,
        show_cell=True,
        spin=spin,
    )

    showmol(viewer, height=VIEWER_HEIGHT, width=VIEWER_WIDTH)

    sol_type_col, cutoff_col = st.columns(2)
    with sol_type_col:
        sol_type = st.selectbox(
            "Solvent Type",
            ["Free Only", "All"],
        )
    with cutoff_col:
        cutoff = st.number_input("threshold", min_value=0.0, max_value=1.0, value=0.5, step=0.1, key="cutoff")

    pred = st.button("🏃‍♂️ Start Cleaning Your MOF")

    upload_id = hashlib.sha256(cif_data).hexdigest()
    checkpoint = "free" if sol_type == "Free Only" else "all"

    if pred:

        safe_name = Path(uploaded_file.name).name
        input_path = input_dir / safe_name
        input_path.write_bytes(cif_data)

        try:
            with st.spinner(f"Running the {checkpoint} solvent model..."):
                runner = get_runner(checkpoint)
                result = runner.clean(
                    cif=input_path,
                    output=output_dir,
                    threshold=cutoff,
                )
        except Exception as error:
            st.error("Your structure without any solvent.")
        else:
            st.session_state["prediction_result"] = result
            st.session_state["prediction_upload_id"] = upload_id
            st.session_state["prediction_checkpoint"] = checkpoint

    result = st.session_state.get("prediction_result")
    result_matches_input = (
        result is not None
        and st.session_state.get("prediction_upload_id") == upload_id
        and st.session_state.get("prediction_checkpoint") == checkpoint
    )

    if result_matches_input:
        colored_header(
            label="Cleaned Structure",
            description="",
            color_name="orange-70",
        )

        available_content = ["Framework"]
        if result.get("solvent") and Path(result["solvent"]).is_file():
            available_content.append("Solvent")
        if st.session_state.get("clean_co") not in available_content:
            st.session_state["clean_co"] = "Framework"

        content_col, type_col, label_col, spin_col = st.columns(4)
        with content_col:
            cleaned_content = st.selectbox(
                "Content",
                available_content,
                key="clean_co",
            )

        with type_col:
            cleaned_render_style = st.selectbox(
                "Style",
                ["Ball and stick", "Stick", "Space filling", "Line"],
                key="clean_rs",
            )
        with label_col:
            cleaned_show_labels = st.selectbox(
                "Element",
                [False, True],
                format_func=lambda value: "Show" if value else "Hide",
                key="clean_sl",
            )
        with spin_col:
            cleaned_spin = st.selectbox(
                "Rotation",
                [False, True],
                format_func=lambda value: "On" if value else "Off",
                key="clean_sp",
            )

        framework_path = Path(result["framework"])
        st.download_button(
            label="⬇️ Download Cleaned MOF",
            data=framework_path.read_bytes(),
            file_name=f"{Path(uploaded_file.name).stem}_gtsr.cif",
            mime="chemical/x-cif",
            on_click="ignore",
        )

        if cleaned_content == "Framework":
            show_cleaned_structure(
                result["framework"],
                cleaned_render_style,
                cleaned_show_labels,
                cleaned_spin,
                "strucinfo_frame",
            )
        elif cleaned_content == "Solvent":
            show_cleaned_structure(
                result["solvent"],
                cleaned_render_style,
                cleaned_show_labels,
                cleaned_spin,
                "strucinfo_sol",
            )
            solvent_smiles = result["solvent_smiles"]
            for i, sol_smi in enumerate(solvent_smiles):
                st.info("Solvents "+ str(i+1) + ": " + sol_smi)

        colored_header(
            label="Activation Stability Prediction",
            description="",
            color_name="orange-70",
        )

        stability_pred = st.button("🏃‍♂️ Want to check stability?")

        if stability_pred:
            with st.spinner("Predicting activation stability..."):
                sta_runner = get_runner(checkpoint="stability")
                score = sta_runner.stability(Path(result["framework"]))

            if score == 1:
                st.info("The cleaned structure is stable.")
            elif score == 0:
                st.error("The cleaned structure is not stable.")
            

    elif result is not None:
        st.info("The target was changed, please re-predict for it.")


colored_header(
    label="Citation",
    description="",
    color_name="orange-70",
)

st.markdown(
    """
```bibtex
@article{gtsr-xyl-group,
  title   = {GTSR: A GNN Based Tool for Solvent Removal from MOF with stability check},
  author  = {Liang, Kairui and Zhao, Guobin and Li, Xiao-Yan},
  journal = {},
  year    = {2026},
  volume  = {},
  number  = {},
  pages   = {},
  doi     = {}
}
```
"""
)
