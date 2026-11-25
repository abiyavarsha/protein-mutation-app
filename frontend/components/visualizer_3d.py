import os
from typing import Optional
import streamlit.components.v1 as components


def render_3d_viewer(
    ref_pdb_content: Optional[str] = None,
    mut_pdb_content: Optional[str] = None,
    superimposed_mut_content: Optional[str] = None,
    view_mode: str = "superimposed", # "reference", "mutant", "superimposed", "side_by_side"
    style: str = "cartoon", # "cartoon", "stick", "sphere"
    color_scheme: str = "protein", # "protein", "plddt"
    mutation_res_num: Optional[int] = None,
    height: int = 550,
    width: str = "100%"
):
    """
    Renders an interactive 3D protein structure viewer using 3Dmol.js WebGL.
    Supports reference, mutant, and superimposed structures with mutation site highlighting
    and pLDDT confidence coloring.
    """
    # Clean PDB content strings for JavaScript injection
    def clean_js_str(s):
        if not s:
            return ""
        return s.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")

    ref_js = clean_js_str(ref_pdb_content)
    mut_js = clean_js_str(mut_pdb_content)
    sup_js = clean_js_str(superimposed_mut_content)

    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <script src="https://3dmol.org/build/3Dmol-min.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                overflow: hidden;
                background-color: #0F172A;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            }}
            #container {{
                width: 100%;
                height: {height}px;
                position: relative;
            }}
            .overlay-badge {{
                position: absolute;
                top: 10px;
                left: 10px;
                background: rgba(15, 23, 42, 0.85);
                backdrop-filter: blur(4px);
                border: 1px solid #334155;
                border-radius: 6px;
                padding: 8px 14px;
                color: #F8FAFC;
                font-size: 12px;
                z-index: 100;
                line-height: 1.4;
                pointer-events: none;
            }}
            .legend-item {{
                display: flex;
                align-items: center;
                gap: 6px;
                margin-top: 2px;
            }}
            .legend-dot {{
                width: 10px;
                height: 10px;
                border-radius: 50%;
                display: inline-block;
            }}
            .controls-hint {{
                position: absolute;
                bottom: 10px;
                right: 10px;
                background: rgba(15, 23, 42, 0.75);
                border: 1px solid #334155;
                border-radius: 4px;
                padding: 4px 8px;
                color: #94A3B8;
                font-size: 11px;
                z-index: 100;
            }}
        </style>
    </head>
    <body>
        <div id="container">
            <div class="overlay-badge">
                <div style="font-weight: 600; margin-bottom: 4px;">3D Structural View ({view_mode.replace('_', ' ').title()})</div>
                <div class="legend-item"><span class="legend-dot" style="background: #00D2D3;"></span> Reference Protein (Cyan)</div>
                <div class="legend-item"><span class="legend-dot" style="background: #FF4757;"></span> Mutant Protein (Magenta / Coral)</div>
                {f'<div class="legend-item"><span class="legend-dot" style="background: #FED330;"></span> Mutation Site (Residue {mutation_res_num})</div>' if mutation_res_num else ''}
            </div>
            <div class="controls-hint">Rotate: Left-Drag | Zoom: Scroll / Right-Drag | Pan: Middle-Drag</div>
        </div>

        <script>
            document.addEventListener("DOMContentLoaded", function() {{
                const element = document.getElementById("container");
                const config = {{ backgroundColor: "#0F172A" }};
                const viewer = $3Dmol.createViewer(element, config);

                const viewMode = "{view_mode}";
                const styleType = "{style}";
                const colorScheme = "{color_scheme}";
                const mutRes = {mutation_res_num if mutation_res_num is not None else "null"};

                const refData = `{ref_js}`;
                const mutData = `{mut_js}`;
                const supData = `{sup_js}`;

                // Helper for pLDDT gradient coloring
                const plddtGrad = function(atom) {{
                    const b = atom.b;
                    if (b >= 90) return "#00539C"; // Very high: Dark blue
                    if (b >= 70) return "#65AFFF"; // Confident: Light blue
                    if (b >= 50) return "#FFDF00"; // Low: Yellow
                    return "#FF4500";             // Very low: Orange/Red
                }};

                if (viewMode === "reference" && refData) {{
                    viewer.addModel(refData, "pdb");
                    if (colorScheme === "plddt") {{
                        viewer.setStyle({{}}, {{ [styleType]: {{ colorfunc: plddtGrad }} }});
                    }} else {{
                        viewer.setStyle({{}}, {{ [styleType]: {{ color: "#00D2D3" }} }});
                    }}
                }} else if (viewMode === "mutant" && mutData) {{
                    viewer.addModel(mutData, "pdb");
                    if (colorScheme === "plddt") {{
                        viewer.setStyle({{}}, {{ [styleType]: {{ colorfunc: plddtGrad }} }});
                    }} else {{
                        viewer.setStyle({{}}, {{ [styleType]: {{ color: "#FF4757" }} }});
                    }}
                }} else if (viewMode === "superimposed") {{
                    // Model 0: Reference
                    if (refData) {{
                        const mRef = viewer.addModel(refData, "pdb");
                        viewer.setStyle({{ model: mRef }}, {{ [styleType]: {{ color: "#00D2D3", opacity: 0.85 }} }});
                    }}
                    // Model 1: Superimposed Mutant
                    const activeMutData = supData || mutData;
                    if (activeMutData) {{
                        const mMut = viewer.addModel(activeMutData, "pdb");
                        viewer.setStyle({{ model: mMut }}, {{ [styleType]: {{ color: "#FF4757", opacity: 0.95 }} }});
                    }}
                }}

                // Highlight mutation residue if provided
                if (mutRes !== null) {{
                    viewer.addStyle({{ resi: mutRes }}, {{
                        sphere: {{ color: "#FED330", radius: 1.2, opacity: 0.9 }},
                        stick: {{ color: "#FED330", radius: 0.35 }}
                    }});
                    viewer.addLabel("Mut: " + mutRes, {{
                        position: {{ resi: mutRes }},
                        backgroundColor: "#1E293B",
                        fontColor: "#FED330",
                        fontSize: 12,
                        backgroundOpacity: 0.8
                    }});
                }}

                viewer.zoomTo();
                viewer.render();
                viewer.spin(true);
                setTimeout(function() {{ viewer.spin(false); }}, 2000);
            }});
        </script>
    </body>
    </html>
    """
    components.html(html_code, height=height, width=None)
