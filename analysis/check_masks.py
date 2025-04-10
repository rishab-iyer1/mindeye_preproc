import argparse
import glob
import hashlib
import nibabel as nib
import numpy as np
import os

def md5_array(arr):
    """Compute md5 hash of a numpy array's raw data."""
    return hashlib.md5(arr.tobytes()).hexdigest()

def main(subject_id, save_final_mask=False, dry_run=False):
    # Find all *_nsdgeneral.nii.gz files
    nii_files = sorted(glob.glob("*_nsdgeneral.nii.gz"))
    
    if len(nii_files) < 2:
        print("❗ Need at least two *_nsdgeneral.nii.gz files to compare.")
        return

    print(f"🔍 Found {len(nii_files)} matching nsdgeneral files:")
    for f in nii_files:
        print(f"  - {f}")

    # Load reference file
    ref_img = nib.load(nii_files[0])
    ref_data = ref_img.get_fdata()
    ref_shape = ref_data.shape
    ref_affine = ref_img.affine
    ref_md5 = md5_array(ref_data)

    print(f"\n🧠 Reference file: {nii_files[0]}")
    print(f"  Shape: {ref_shape}")

    # Check all other files
    for f in nii_files[1:]:
        img = nib.load(f)
        data = img.get_fdata()
        shape = data.shape
        affine = img.affine
        md5 = md5_array(data)

        if shape != ref_shape:
            raise ValueError(f"❌ Shape mismatch in {f}: {shape} vs {ref_shape}")
        if not np.allclose(affine, ref_affine):
            raise ValueError(f"❌ Affine mismatch in {f}")
        if md5 != ref_md5:
            raise ValueError(f"❌ Data mismatch in {f}")
        
        print(f"  ✔ {f}: shape, affine, and data match")
    
    outname = os.path.abspath(f"sub-{subject_id}_final_nsdgeneral.nii.gz")
    # Dry run diagnostics
    if dry_run:
        print("\n🔎 Dry run diagnostics:")
        for i, f in enumerate(nii_files):
            img = nib.load(f)
            data = img.get_fdata()
            md5 = md5_array(data)
            nonzero = int(np.count_nonzero(data))
            print(f"  File {i+1}: {f}")
            print(f"    - MD5: {md5}")
            print(f"    - Nonzero voxels: {nonzero:,}")
        
        print(f"  Save path: {outname}")
        print("🛑 Dry run enabled — no output saved.")
        return

    # Save the final mask if requested
    if save_final_mask:
        nib.save(ref_img, outname)
        print(f"\n✅ All files match. Final nsdgeneral mask saved as: {outname}")
    else:
        print(f"  Save path: {outname}")
        print("\n✅ All files match. No output saved (use --save-final-mask to save result).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify and optionally save *_nsdgeneral.nii.gz consistency")
    parser.add_argument("subject_id", type=str, help="Subject ID (e.g., 01)")
    parser.add_argument("--save-final-mask", action="store_true", help="Save final verified nsdgeneral mask")
    parser.add_argument("--dry-run", action="store_true", help="Print diagnostics but do not save anything")
    args = parser.parse_args()

    main(subject_id=args.subject_id, save_final_mask=args.save_final_mask, dry_run=args.dry_run)
