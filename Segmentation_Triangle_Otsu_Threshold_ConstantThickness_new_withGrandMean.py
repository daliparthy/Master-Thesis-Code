# python q:\code\process_sem_final4.py

import cv2
import os
import glob
import numpy as np
import pandas as pd
import re
# =============================================================================
# CONFIGURATION
# =============================================================================

input_dir = r"Q:\Gemini SEM\Master Arbeit\5nm\Ar\PlainImages\seg5k\Processed_Images"

pixel_size_nm = 22.33

# Minimum black hole/substrate component size to KEEP.
# Black regions smaller than this will be removed/filled.
min_hole_size = 1

# Optional: minimum white particle size to KEEP.
# Set to 1 or None if you do not want to remove small white particles.
min_particle_size = 2

output_dir = os.path.join(input_dir, "Processed_Images3")
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def pct_white(mask):
    """Return percentage of white pixels in a binary image."""
    return 100 * cv2.countNonZero(mask) / mask.size

""" Use for circularish islands / dark holes on a fairly uniform grey film """
def apply_user_threshold_logic(img, filename="Unknown", min_hole_size=30, min_hole_darkness=90, min_particle_size=1, debug=False):

    # Global threshold (Triangle, with Otsu fallback) to prevent 'white-out' images.
    # This works well when film vs. hole brightness is genuinely bimodal
    # (a single grey film level + a darker hole level), which is the case
    # for these SEM images.

    # Convert to grayscale if needed
    if len(img.shape) == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Convert to 8-bit if needed
    if img.dtype != np.uint8:
        img = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)

    # --- Apply CLAHE to normalize contrast across the batch ---
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    img_normalized = clahe.apply(img)

    # Denoising: median blur removes sensor grain without smearing edges
    blurred = cv2.medianBlur(img_normalized, 3)

    # 1. Try Triangle Threshold first (good for a small, separated dark mode)
    ret, thresh = cv2.threshold(
        blurred,
        0,
        255,
        cv2.THRESH_BINARY + cv2.THRESH_TRIANGLE
    )

    # 2. Sanity check: Triangle can fail on some histogram shapes.
    #    Fall back to Otsu if the result looks degenerate (either the
    #    threshold value itself is implausibly low, or the resulting mask
    #    is almost entirely black/white, which signals a bad split).
    white_pct = pct_white(thresh)
    if ret < 10 or white_pct < 2 or white_pct > 98:
        ret_otsu, thresh_otsu = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        if debug:
            print(f"--- Fallback Triggered for {filename} ---")
            print(f"Triangle result rejected (thresh={ret}, white%={white_pct:.1f}); using Otsu ({ret_otsu})")

        ret = ret_otsu
        thresh = thresh_otsu

    binary_initial = thresh.copy()

    # STAGE A: REMOVE SMALL BLACK HOLES
    inv_black_regions = cv2.bitwise_not(binary_initial)
    n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(inv_black_regions, connectivity=4)

    valid_black_labels = np.where(stats[:, cv2.CC_STAT_AREA] >= min_hole_size)[0]
    valid_black_labels = valid_black_labels[valid_black_labels > 0]

    cleaned_inv_black_regions = np.isin(labels, valid_black_labels).astype(np.uint8) * 255
    binary_clean_holes = cv2.bitwise_not(cleaned_inv_black_regions)

    # STAGE B: OPTIONAL REMOVE SMALL WHITE PARTICLES
    if min_particle_size is not None and min_particle_size > 1:
        n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(
            binary_clean_holes,
            connectivity=8
        )
        valid_particle_labels = np.where(stats[:, cv2.CC_STAT_AREA] >= min_particle_size)[0]
        valid_particle_labels = valid_particle_labels[valid_particle_labels > 0]
        final_binary = np.isin(labels, valid_particle_labels).astype(np.uint8) * 255
    else:
        final_binary = binary_clean_holes

    if debug:
        print(f"---- threshold debug for {filename} ----")
        print(f"Final threshold value: {ret}")
        print(f"final_binary white %: {pct_white(final_binary):.2f}")
        print("-------------------------")

    return final_binary, ret, binary_initial, binary_clean_holes


def is_on_edge(cnt, shape):
    h, w = shape[:2]
    pts = cnt[:, 0, :]

    return (
        np.any(pts[:, 0] <= 1) or
        np.any(pts[:, 0] >= w - 2) or
        np.any(pts[:, 1] <= 1) or
        np.any(pts[:, 1] >= h - 2)
    )


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def process_images(input_dir, output_dir, px_nm):
    image_paths = (
        glob.glob(os.path.join(input_dir, "*.png")) +
        glob.glob(os.path.join(input_dir, "*.tiff"))
    )

    if not image_paths:
        print("No images found.")
        return

    all_results = []
    
    # --- GLOBAL LISTS FOR GRAND STATISTICS ---
    global_diams = []
    global_circs = []
    global_densities = []
    global_ca=[]

    for img_path in image_paths:
        filename = os.path.basename(img_path)

        img = cv2.imread(img_path, cv2.IMREAD_UNCHANGED)
        if img is None:
            print(f"Could not read {filename}, skipping.")
            continue

        # ---------------------------------------------------------------------
        # 1. Generate binary masks
        # final_binary: white = particles/film, black = holes/substrate/background
        # binary_initial: raw threshold before removing small black holes
        # binary_clean_holes: after black-hole cleanup, before white-particle cleanup
        # ---------------------------------------------------------------------

        binary, thresh_val, binary_initial, binary_clean_holes = apply_user_threshold_logic(
            img, filename=filename,
            min_hole_size=min_hole_size,
            min_particle_size=min_particle_size,
            debug=True
        )

        # ---------------------------------------------------------------------
        # 2. Image metrics
        # ---------------------------------------------------------------------
        total_px = binary.shape[0] * binary.shape[1]

        island_px_all = cv2.countNonZero(binary)
        pixel_coverage_pct = (island_px_all / total_px) * 100

        pixel_size_um = px_nm * 1e-3
        image_area_um2 = (
            binary.shape[0] * pixel_size_um
        ) * (
            binary.shape[1] * pixel_size_um
        )

        # ---------------------------------------------------------------------
        # 3. Contour detection
        # ---------------------------------------------------------------------
        contours, hierarchy = cv2.findContours(
            binary,
            cv2.RETR_CCOMP,
            cv2.CHAIN_APPROX_SIMPLE
        )

        img_8u = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)

        if len(img_8u.shape) == 2:
            overlay = cv2.cvtColor(img_8u, cv2.COLOR_GRAY2BGR)
        else:
            overlay = img_8u.copy()

        thickness = max(1, round(min(binary.shape[:2]) / 512))

        diams_nm = []
        circularities = []
        volumes_nm3 = []
        surface_areas_nm2 = []

        total_contour_length_nm = 0.0
        total_interior_area_px = 0.0

        if hierarchy is not None:
            for i, cnt in enumerate(contours):
                if hierarchy[0, i, 3] == -1:
                    if is_on_edge(cnt, binary.shape):
                        continue

                    net_area_px = cv2.contourArea(cnt)
                    perimeter = cv2.arcLength(cnt, True)

                    child_idx = hierarchy[0, i, 2]
                    while child_idx != -1:
                        net_area_px -= cv2.contourArea(contours[child_idx])
                        cv2.drawContours(overlay, contours, child_idx, (0, 0, 255), thickness)
                        child_idx = hierarchy[0, child_idx][0]

                    if net_area_px > 0:
                        d_nm = np.sqrt((4 * net_area_px) / np.pi) * px_nm
                        circ = (4 * np.pi * net_area_px) / (perimeter ** 2) if perimeter > 0 else 0

                        diams_nm.append(d_nm)
                        circularities.append(circ)
                        volumes_nm3.append((np.pi / 6) * (d_nm ** 3))
                        surface_areas_nm2.append(np.pi * (d_nm ** 2))

                        total_contour_length_nm += perimeter * px_nm
                        total_interior_area_px += net_area_px

                        cv2.drawContours(overlay, contours, i, (0, 0, 255), thickness)
        
        # --- APPEND TO GLOBAL LISTS ---
        global_diams.extend(diams_nm)
        global_circs.extend(circularities)

        # ---------------------------------------------------------------------
        # 4. Statistics
        # ---------------------------------------------------------------------
        count = len(diams_nm)

        if count > 0:
            mean_d = np.mean(diams_nm)
            std_d = np.std(diams_nm)
            mean_c = np.mean(circularities)
            std_c = np.std(circularities)
            tot_vol = np.sum(volumes_nm3)
            tot_surf = np.sum(surface_areas_nm2)
            s_v_ratio = tot_surf / tot_vol if tot_vol > 0 else 0
        else:
            mean_d = 0
            std_d = 0
            mean_c = 0
            std_c = 0
            tot_vol = 0
            tot_surf = 0
            s_v_ratio = 0

        interior_coverage_pct = (total_interior_area_px / total_px) * 100
        particle_density = count / image_area_um2 if image_area_um2 > 0 else 0
        norm_contour_length = total_contour_length_nm / image_area_um2 if image_area_um2 > 0 else 0

        global_densities.append(particle_density)
        global_ca.append(pixel_coverage_pct)
        
        all_results.append([
            filename,
            mean_d,
            std_d,
            mean_c,
            std_c,
            count,
            particle_density,
            total_contour_length_nm,
            norm_contour_length,
            pixel_coverage_pct,
            interior_coverage_pct,
            tot_vol,
            tot_surf,
            s_v_ratio
        ])

        # ---------------------------------------------------------------------
        # 5. Save visual outputs
        # ---------------------------------------------------------------------
        clean_filename = os.path.splitext(filename)[0]

        binary_path = os.path.join(output_dir, f"binary_{clean_filename}.png")
        cv2.imwrite(binary_path, binary)

        overlay_path = os.path.join(output_dir, f"overlay_{clean_filename}.png")
        success = cv2.imwrite(overlay_path, overlay)

        if not success:
            print(f"FAILED to save overlay: {overlay_path}")
        else:
            print(f"Processed: {clean_filename} | Islands Analyzed: {count}")

    # -------------------------------------------------------------------------
    # 6. Calculate Grand Statistics & Export Excel results
    # -------------------------------------------------------------------------
    cols = [
        "Image",
        "Mean Diameter (nm)",
        "Std Dev Diameter (nm)",
        "Mean Circularity",
        "Std Dev Circularity",
        "Particle Count",
        "Particle Density (particles/µm²)",
        "Total Contour Length (nm)",
        "Total Contour Length Normalised (nm/µm²)",
        "Covered Area – all pixels incl. edge (%)",
        "Covered Area – interior particles only (%)",
        "Total Volume (nm³)",
        "Total Surface Area (nm²)",
        "Surface-to-Volume Ratio (nm⁻¹)"
    ]

    df = pd.DataFrame(all_results, columns=cols)

    # --- GRAND STATISTICS CALCULATIONS ---
    grand_mean_ecd = np.mean(global_diams) if global_diams else 0
    grand_std_ecd = np.std(global_diams) if global_diams else 0
    
    grand_mean_circ = np.mean(global_circs) if global_circs else 0
    grand_std_circ = np.std(global_circs) if global_circs else 0
    
    grand_mean_density = np.mean(global_densities) if global_densities else 0
    grand_std_density = np.std(global_densities) if global_densities else 0

    grand_mean_ca = np.mean(global_ca) if global_ca else 0
    grand_std_ca = np.std(global_ca) if global_ca else 0

    print("\n" + "="*50)
    print("GRAND STATISTICS (Across all processed images)")
    print("="*50)
    print(f"Grand Mean ECD: {grand_mean_ecd:.2f} nm ± {grand_std_ecd:.2f} nm")
    print(f"Grand Mean Circularity: {grand_mean_circ:.4f} ± {grand_std_circ:.4f}")
    print(f"Grand Mean Particle Density: {grand_mean_density:.2f} ± {grand_std_density:.2f} (particles/µm²)")
    print(f"Grand Mean Covered Area: {grand_mean_ca:.2f} ± {grand_std_ca:.2f}")
    print("="*50 + "\n")

    # Append Grand Stats to the DataFrame
    empty_row = {col: None for col in cols}
    
    mean_row = empty_row.copy()
    mean_row["Image"] = "GRAND MEAN"
    mean_row["Mean Diameter (nm)"] = grand_mean_ecd
    mean_row["Mean Circularity"] = grand_mean_circ
    mean_row["Particle Density (particles/µm²)"] = grand_mean_density
    mean_row["Covered Area – all pixels incl. edge (%)"] = grand_mean_ca

    std_row = empty_row.copy()
    std_row["Image"] = "GRAND STD DEV"
    std_row["Mean Diameter (nm)"] = grand_std_ecd
    std_row["Mean Circularity"] = grand_std_circ
    std_row["Particle Density (particles/µm²)"] = grand_std_density
    std_row["Covered Area – all pixels incl. edge (%)"] = grand_std_ca

    summary_df = pd.DataFrame([mean_row, std_row])
    df = pd.concat([df, summary_df], ignore_index=True)

    excel_path = os.path.join(output_dir, "SEM_Advanced_Analysis.xlsx")
    df.to_excel(excel_path, index=False)

    print(f"Saved Excel results (with Grand Stats appended) to: {excel_path}")


if __name__ == "__main__":
    process_images(
        input_dir=input_dir,
        output_dir=output_dir,
        px_nm=pixel_size_nm,
    )