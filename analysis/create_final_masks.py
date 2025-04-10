import argparse
import glob
import os
import nibabel as nib
import numpy as np
from nilearn.masking import intersect_masks

def main(subject_id, sessions, dry_run=True):
    # Find all brain masks for the specified sessions
    all_masks = []
    for sess in sessions:
        matches = glob.glob(f"*_{sess}_*_brain.nii.gz")
        filtered = [f for f in matches if f"_ses-" in f or f"_{sess}_" in f]
        all_masks.extend(filtered)

    if not all_masks:
        print("No matching brain masks found.")
        return

    print(f"Found {len(all_masks)} mask(s):")
    for f in all_masks:
        print(f"  {f}")

    # Load masks and check shape/affine consistency
    reference_img = nib.load(all_masks[0])
    ref_shape = reference_img.shape
    ref_affine = reference_img.affine
    print(f"Reference shape: {ref_shape}")
    print("Checking shape and affine consistency...")

    mask_imgs = [reference_img]

    for f in all_masks[1:]:
        img = nib.load(f)
        if img.shape != ref_shape:
            raise ValueError(f"Shape mismatch in {f}: {img.shape} vs {ref_shape}")
        if not np.allclose(img.affine, ref_affine):
            raise ValueError(f"Affine mismatch in {f}")
        print(f"  ✔ {f}: shape and affine match")
        mask_imgs.append(img)

    # Intersect masks (most conservative = only voxels present in all masks)
    intersected = intersect_masks(mask_imgs, threshold=1.0, connected=False)

    # Binarize result
    binarized_data = (intersected.get_fdata() > 0).astype(np.uint8)
    bin_img = nib.Nifti1Image(binarized_data, affine=ref_affine)

    output_name = os.path.abspath(f"sub-{subject_id}_final_brain.nii.gz")
    if dry_run:
        print("\n🔍 Dry run diagnostics:")
        for i, img in enumerate(mask_imgs):
            nvox = int(np.count_nonzero(img.get_fdata()))
            print(f"  Mask {i+1}: {nvox:,} nonzero voxels")
        final_voxels = int(np.count_nonzero(binarized_data))
        pct = (final_voxels / np.count_nonzero(mask_imgs[0].get_fdata())) * 100
        print(f"  Intersected mask: {final_voxels:,} nonzero voxels ({pct:.2f}% of first mask)")
        print(f"  Save path: {output_name}")
        print("🛑 Dry run enabled — final mask NOT saved.")
    else:
        nib.save(bin_img, output_name)
        print(f"✅ Saved final brain mask as: {output_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Intersect brain masks conservatively")
    parser.add_argument("subject_id", type=str, help="Subject ID (e.g., 01)")
    parser.add_argument("sessions", nargs="+", help="List of session names (e.g., ses-01 ses-02)")
    parser.add_argument("--dry-run", action="store_true", help="Run without saving the final mask")


    args = parser.parse_args()
    main(args.subject_id, args.sessions, dry_run=args.dry_run)

