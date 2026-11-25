import os
import time
import urllib.request
import urllib.error
from typing import Dict, Any, Optional, Tuple
from Bio.PDB import PDBParser


ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_ESMFOLD_SEQ_LENGTH = 400


class StructureService:

    @staticmethod
    def validate_pdb(pdb_content_or_path: str) -> Tuple[bool, str, int]:
        """
        Validates that a PDB string or file contains valid PDB format with ATOM records.
        Returns: (is_valid, message, atom_count)
        """
        lines = []
        if os.path.exists(pdb_content_or_path):
            with open(pdb_content_or_path, "r") as f:
                lines = f.readlines()
        else:
            lines = pdb_content_or_path.splitlines()
            
        atom_count = 0
        ca_count = 0
        for line in lines:
            if line.startswith("ATOM  "):
                atom_count += 1
                if len(line) >= 16 and line[12:16].strip() == "CA":
                    ca_count += 1
                    
        if atom_count == 0 or ca_count == 0:
            return False, f"Invalid PDB: Found {atom_count} ATOM records and {ca_count} C-alpha atoms.", 0
            
        return True, f"Valid PDB: {atom_count} atoms, {ca_count} C-alpha residues.", ca_count

    @staticmethod
    def predict_with_esmfold(
        sequence: str,
        output_path: str,
        max_retries: int = 3,
        timeout: int = 90
    ) -> Dict[str, Any]:
        """
        Submits protein sequence to the hosted ESMFold API and saves the validated PDB.
        Handles network errors, rate limits, timeouts, and HTTP status codes.
        Never fabricates a structure on failure.
        """
        clean_seq = sequence.strip().upper()
        
        if len(clean_seq) > MAX_ESMFOLD_SEQ_LENGTH:
            raise ValueError(
                f"Sequence length ({len(clean_seq)} aa) exceeds hosted ESMFold API limit "
                f"of {MAX_ESMFOLD_SEQ_LENGTH} aa (HTTP 413: Request Entity Too Large). "
                f"Please select a functional domain or window containing the mutation site."
            )
            
        t0 = time.time()
        last_error = None
        
        for attempt in range(1, max_retries + 1):
            try:
                req = urllib.request.Request(
                    ESMFOLD_API_URL,
                    data=clean_seq.encode("utf-8"),
                    headers={
                        "Content-Type": "text/plain",
                        "User-Agent": "ProteinMutationAnalysisApp/1.0"
                    }
                )
                
                with urllib.request.urlopen(req, timeout=timeout) as response:
                    if response.status == 200:
                        pdb_content = response.read().decode("utf-8")
                        
                        # Validate returned content
                        is_valid, msg, ca_count = StructureService.validate_pdb(pdb_content)
                        if not is_valid:
                            raise ValueError(f"ESMFold API returned non-PDB content: {msg}")
                            
                        # Save structure to filesystem
                        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
                        with open(output_path, "w") as f:
                            f.write(pdb_content)
                            
                        elapsed = time.time() - t0
                        return {
                            "success": True,
                            "service": "ESMFold",
                            "output_path": output_path,
                            "sequence_length": len(clean_seq),
                            "ca_atom_count": ca_count,
                            "elapsed_seconds": round(elapsed, 2),
                            "attempt": attempt,
                            "message": f"Structure successfully predicted by ESMFold in {elapsed:.1f}s."
                        }
                    else:
                        raise RuntimeError(f"ESMFold API returned HTTP status {response.status}")
                        
            except urllib.error.HTTPError as he:
                if he.code == 504:
                    last_error = "HTTP 504: Gateway Timeout (Hosted ESMFold public server congested)"
                elif he.code == 429:
                    last_error = "HTTP 429: Too Many Requests (ESMFold rate limit reached)"
                elif he.code == 413:
                    last_error = "HTTP 413: Sequence Too Large (> 400 aa)"
                    break
                else:
                    last_error = f"HTTP {he.code}: {he.reason}"
                time.sleep(1.5 * attempt)
            except Exception as e:
                last_error = str(e)
                time.sleep(1.5 * attempt)
                
        # Failure
        raise RuntimeError(
            f"Hosted ESMFold API service temporarily unavailable ({last_error}). "
            f"Public server (api.esmatlas.com) may be congested. Please try again or provide a local PDB."
        )

    @staticmethod
    def get_structure_for_protein(
        sequence: str,
        protein_name: str,
        structure_type: str, # "reference" or mutant label, e.g. "E6V_mutant"
        base_dir: str,
        force_predict: bool = False
    ) -> Dict[str, Any]:
        """
        Retrieves or generates tertiary structure for a given sequence and protein.
        Checks local filesystem cache in data/structures/{protein}/{structure_type}.pdb.
        If not cached, requests prediction from ESMFold (or domain window if > 400 aa).
        """
        struct_dir = os.path.join(base_dir, "data", "structures", protein_name)
        os.makedirs(struct_dir, exist_ok=True)
        
        target_filename = f"{structure_type}.pdb" if not structure_type.endswith(".pdb") else structure_type
        target_path = os.path.join(struct_dir, target_filename)
        
        # Check cache
        if not force_predict and os.path.exists(target_path):
            is_valid, msg, ca_count = StructureService.validate_pdb(target_path)
            if is_valid:
                return {
                    "success": True,
                    "service": "Local Cache",
                    "output_path": target_path,
                    "sequence_length": len(sequence),
                    "ca_atom_count": ca_count,
                    "elapsed_seconds": 0.0,
                    "message": f"Loaded validated structure from local cache: {target_filename}"
                }
                
        # If sequence is <= 400 aa, run ESMFold directly
        if len(sequence) <= MAX_ESMFOLD_SEQ_LENGTH:
            return StructureService.predict_with_esmfold(sequence, target_path)
        else:
            # For sequences > 400 aa, determine appropriate domain or notify
            raise ValueError(
                f"Protein {protein_name} sequence length ({len(sequence)} aa) exceeds hosted "
                f"ESMFold API limit of {MAX_ESMFOLD_SEQ_LENGTH} aa. "
                f"Please fold a defined domain (e.g. NBD1 for CFTR, Kinase domain for EGFR)."
            )
