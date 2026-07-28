Graph A represents a **"Methane Tug-of-War"** balance. To create it, we are using two specific tabs from your Excel workbook:

1. **`pmoA_OTUs_rel_abundance`** (Methane Eaters)
2. **`mcrA_OTUs_rel_abundance`** (Methane Producers)

---

### Step 1: The Rows (How the Data Matches Up)

Every single row in your spreadsheet represents a **unique biological sample** taken from a specific mud pot in the experiment.

If you look at the very first column, labeled **`Group`**, you will notice that both sheets have the exact same sample IDs in the exact same order.

For example, look at row #4 in both sheets:

* In the $pmoA$ sheet, row 4 is `IR64_Control_R1` (at the `Early_tillering_stage`).
* In the $mcrA$ sheet, row 4 is also `IR64_Control_R1` (at the `Early_tillering_stage`).

Because the rows match perfectly, the code reads row 4 from both sheets simultaneously to calculate a balance for that exact plant plot.

---

### Step 2: The Columns (What is Being Calculated)

#### 1. The Metadata Columns (X-Axis and Split Panels)

The code ignores the first two columns (`Group` and `Stage`) during the math, but uses them to organize the final layout of the graph:

* **`Stage` column:** Tells the graph where to plot the point horizontally along the timeline (**X-Axis**).
* **`Group` column:** The code looks at a name like `IR64_Control_R1` and splits it up to know that the variety is `IR64` and the treatment is `Control`. This dictates the color of the lines and splits the graph into side-by-side comparison panels (`IR64` panel vs. `Nipponbare` panel).

#### 2. The Abundance Columns (The Math)

Every column after the first two is an OTU (a specific microbial species).

* In the $pmoA$ sheet, you have columns running from `Otu001` all the way to `Otu122`.
* In the $mcrA$ sheet, you have columns running from `Otu01` to `Otu18`.

The values inside these columns are decimal percentages (Relative Abundance). For example, if a cell says `0.50`, it means that specific microbe makes up 50% of that community.

---

### Step 3: How the Math Creates the Y-Axis Value

To get a single number that tells us who is winning the environmental "Tug-of-War," the code handles the columns like this:

1. **Look at the $pmoA$ row columns:** It looks across all 122 OTU columns for a single sample row and calculates its **Shannon Diversity Index**. This index looks at how many different methane-eating species are present and how evenly distributed their percentages are. Let's say it gets a score of **2.4**.
2. **Look at the $mcrA$ row columns:** It looks across all 18 OTU columns for the exact same sample row and calculates its diversity score. Let's say it gets a score of **1.2**.
3. **Divide them (The Ratio):** The code divides the $pmoA$ score by the $mcrA$ score ($2.4 / 1.2 = 2.0$) and takes the $log_2$ of that number.

This final calculated number becomes the exact spot on the **Y-Axis** for that sample!

### Summary of Graph A Layout:

* **Row data used for:** Matching individual samples across both microbial worlds.
* **OTU columns used for:** Calculating the community complexity scores (Diversity) to see which community is dominant.
* **Stage column used for:** The X-Axis timeline tracking.
* **Group column used for:** Color coding the lines (`Control` vs `KH32C`) and sorting them into separate sub-plots (`IR64` vs `Nipponbare`).

---------------

