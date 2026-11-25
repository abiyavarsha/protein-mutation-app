import os
import datetime
from typing import Dict, Any, List, Optional
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether, HRFlowable
)


class ReportService:

    @staticmethod
    def generate_pdf_report(
        analysis_data: Dict[str, Any],
        output_pdf_path: str
    ) -> str:
        """
        Generates a comprehensive, publication-grade PDF report summarizing:
        - Protein identification and UniProt info
        - Sequence alignment metrics and detected mutations
        - Structure prediction pLDDT confidence
        - Kabsch superimposition metrics (Global RMSD, Local RMSD, Coverage, TM-score)
        - Biophysical property deltas
        - Scientific interpretation with clear disclaimers
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_pdf_path)), exist_ok=True)
        
        doc = SimpleDocTemplate(
            output_pdf_path,
            pagesize=letter,
            rightMargin=40,
            leftMargin=40,
            topMargin=40,
            bottomMargin=40
        )
        
        styles = getSampleStyleSheet()
        
        # Custom styles
        title_style = ParagraphStyle(
            "DocTitle",
            parent=styles["Heading1"],
            fontSize=20,
            leading=24,
            textColor=colors.HexColor("#1A365D"),
            fontName="Helvetica-Bold",
            spaceAfter=6
        )
        subtitle_style = ParagraphStyle(
            "DocSubtitle",
            parent=styles["Normal"],
            fontSize=10,
            leading=14,
            textColor=colors.HexColor("#4A5568"),
            spaceAfter=15
        )
        h2_style = ParagraphStyle(
            "H2",
            parent=styles["Heading2"],
            fontSize=13,
            leading=16,
            textColor=colors.HexColor("#2B6CB0"),
            fontName="Helvetica-Bold",
            spaceBefore=12,
            spaceAfter=6
        )
        body_style = ParagraphStyle(
            "Body",
            parent=styles["Normal"],
            fontSize=9,
            leading=13,
            textColor=colors.HexColor("#2D3748")
        )
        callout_style = ParagraphStyle(
            "Callout",
            parent=styles["Normal"],
            fontSize=8.5,
            leading=12,
            textColor=colors.HexColor("#2C5282")
        )
        disclaimer_style = ParagraphStyle(
            "Disclaimer",
            parent=styles["Normal"],
            fontSize=8,
            leading=11,
            textColor=colors.HexColor("#742A2A")
        )
        
        story = []
        
        # 1. Header Banner
        story.append(Paragraph("Protein Mutation Analysis & Tertiary Structure Comparison", title_style))
        now_str = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        protein_name = analysis_data.get("protein_name", "Unknown")
        mutant_label = analysis_data.get("mutant_label", "Mutant")
        job_id = analysis_data.get("job_id", "N/A")
        
        story.append(Paragraph(
            f"<b>Target Protein:</b> {protein_name} &nbsp;|&nbsp; <b>Mutation:</b> {mutant_label} &nbsp;|&nbsp; "
            f"<b>Job ID:</b> {job_id} &nbsp;|&nbsp; <b>Generated:</b> {now_str}",
            subtitle_style
        ))
        story.append(HRFlowable(width="100%", thickness=1.5, color=colors.HexColor("#2B6CB0"), spaceAfter=12))
        
        # 2. Executive Summary Metrics Table
        story.append(Paragraph("1. Executive Summary & Structural Metrics", h2_style))
        
        global_rmsd = analysis_data.get("global_rmsd", "N/A")
        local_rmsd = analysis_data.get("local_rmsd", "N/A")
        seq_id = analysis_data.get("sequence_identity", "N/A")
        coverage = analysis_data.get("alignment_coverage", "N/A")
        tm_score = analysis_data.get("tm_score", "N/A")
        matched_atoms = analysis_data.get("matched_ca_atoms_count", "N/A")
        
        summary_table_data = [
            [
                Paragraph("<b>Global C-alpha RMSD:</b>", body_style),
                Paragraph(f"<b>{global_rmsd:.4f} Å</b>" if isinstance(global_rmsd, (int, float)) else str(global_rmsd), body_style),
                Paragraph("<b>Local C-alpha RMSD:</b>", body_style),
                Paragraph(f"<b>{local_rmsd:.4f} Å</b>" if isinstance(local_rmsd, (int, float)) else str(local_rmsd), body_style)
            ],
            [
                Paragraph("<b>Sequence Identity:</b>", body_style),
                Paragraph(f"{seq_id:.2f}%" if isinstance(seq_id, (int, float)) else str(seq_id), body_style),
                Paragraph("<b>Structural Coverage:</b>", body_style),
                Paragraph(f"{coverage:.1f}%" if isinstance(coverage, (int, float)) else str(coverage), body_style)
            ],
            [
                Paragraph("<b>Matched CA Atoms:</b>", body_style),
                Paragraph(str(matched_atoms), body_style),
                Paragraph("<b>TM-Score (est.):</b>", body_style),
                Paragraph(f"{tm_score:.4f}" if isinstance(tm_score, (int, float)) else str(tm_score), body_style)
            ]
        ]
        
        t_summary = Table(summary_table_data, colWidths=[130, 130, 130, 130])
        t_summary.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor("#F7FAFC")),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#CBD5E0")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#E2E8F0")),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
        ]))
        story.append(t_summary)
        story.append(Spacer(1, 12))
        
        # 3. Detected Mutations Table
        story.append(Paragraph("2. Sequence Alignment & Mutation Identification", h2_style))
        mutations = analysis_data.get("mutations", [])
        
        if mutations:
            mut_table_data = [
                ["Notation", "Type", "Ref Pos", "Mature Pos", "Ref AA", "Mut AA", "Description"]
            ]
            for m in mutations:
                mut_table_data.append([
                    m.get("canonical_notation", "-"),
                    m.get("mutation_type", "-"),
                    str(m.get("ref_pos", "-")),
                    str(m.get("mature_pos", "-")),
                    m.get("ref_aa", "-") or "-",
                    m.get("mut_aa", "-") or "-",
                    Paragraph(m.get("description", "-"), body_style)
                ])
                
            t_mut = Table(mut_table_data, colWidths=[65, 75, 45, 55, 45, 45, 190])
            t_mut.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#2B6CB0")),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                ('ALIGN', (0, 0), (5, -1), 'CENTER'),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#CBD5E0")),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7FAFC")]),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ]))
            story.append(t_mut)
        else:
            story.append(Paragraph("No mutations identified between sequences.", body_style))
            
        story.append(Spacer(1, 12))
        
        # 4. Biophysical Property Alterations
        story.append(Paragraph("3. Biophysical Property Alterations", h2_style))
        biophysical = analysis_data.get("biophysical_changes", [])
        if biophysical:
            bio_table_data = [
                ["Mutation", "Type", "Ref Residue", "Mut Residue", "ΔHydropathy", "ΔCharge", "Interpretation"]
            ]
            for b in biophysical:
                bio_table_data.append([
                    b.get("mature_notation") or b.get("notation", "-"),
                    b.get("type", "-"),
                    b.get("ref_residue", "-"),
                    b.get("mut_residue", "-"),
                    str(b.get("delta_hydropathy", "N/A")),
                    str(b.get("delta_charge", "N/A")),
                    Paragraph(b.get("interpretation", "-"), body_style)
                ])
            t_bio = Table(bio_table_data, colWidths=[65, 65, 75, 75, 55, 45, 140])
            t_bio.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#4A5568")),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#CBD5E0")),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7FAFC")]),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ]))
            story.append(t_bio)
        story.append(Spacer(1, 12))
        
        # 5. Scientific Interpretation
        story.append(Paragraph("4. Scientific Discussion & Structural Interpretation", h2_style))
        discussion_text = analysis_data.get("discussion", "")
        for line in discussion_text.split("\n"):
            line = line.strip()
            if line:
                # Convert markdown bold to html
                line_html = line.replace("**", "<b>", 1)
                line_html = line_html.replace("**", "</b>", 1)
                story.append(Paragraph(line_html, body_style))
                story.append(Spacer(1, 3))
                
        story.append(Spacer(1, 10))
        
        # 6. Methodological Notes and Limitations
        story.append(Paragraph("5. Methodological Notes & Scientific Disclaimers", h2_style))
        disclaimer_box = [
            [Paragraph(
                "<b>MANDATORY SCIENTIFIC LIMITATIONS:</b><br/>"
                "• <b>RMSD Interpretation:</b> Global and local RMSD metrics quantify geometric root-mean-square deviations between predicted alpha-carbon atom coordinates. "
                "RMSD alone DOES NOT prove protein destabilization, misfolding, or loss of function.<br/>"
                "• <b>Pathogenicity:</b> Pathogenicity cannot be inferred solely from static predicted models. Additional functional assays, clinical cohorts, and population genetics are required.<br/>"
                "• <b>Thermodynamic Stability (ΔΔG):</b> ΔΔG must not be derived from RMSD; independent rigorous thermodynamic methods (e.g. FoldX, Rosetta, or FEP) are required.<br/>"
                "• <b>Prediction Confidence:</b> Tertiary structures predicted by ESMFold v1 provide confidence metrics via per-residue pLDDT scores in the B-factor field. Regions with pLDDT &lt; 70 indicate conformational flexibility or lower prediction confidence.",
                disclaimer_style
            )]
        ]
        t_disc = Table(disclaimer_box, colWidths=[520])
        t_disc.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor("#FFF5F5")),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#FEB2B2")),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('LEFTPADDING', (0, 0), (-1, -1), 8),
            ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ]))
        story.append(t_disc)
        
        # Build PDF
        doc.build(story)
        return output_pdf_path
