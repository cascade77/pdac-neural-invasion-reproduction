## Objectives:
Trying to investigate how a certain gene(TGFBI) certain nerve cells ( Schwann cells) behave /imapcat in pancreatic cancer

-To download TCGA pancreatic cancer patient data and derive wether patients with higher or lower TGFBI gene expression survve longer<br>
-To run cell chat and find out which cells are sending TGF-beta signals to schwannn cells<br>
-To run Monocle3 to find out how schwann cells change their state ver time during  pancreatihow oc cancer<br>

## Background

Pancreatic cancer is one of the most deadly cancers  mainly because it is usually
detected very late and does not respond well to most treatments.

One of the main reasons is neural invasion. This is basically
when cancer cells start to move and invade the nerves running through the pancreas.
Patients whose tumors show high neural invasion show worse results. But befre recently, noone knew exactly which cell types
were responsible for enabling or resisting this nerve invasion.

### Schwann cells:
Wrap around and protect nerves. These are
found along nerves throughout the body. When
cancerous cells invade nerves, Schwann cells are right at the site of invasion.
The paper we are trying to reproduce found that Schwann cells are not just passive bystanders but they
actively interact with cancer cells and other immune cells.

### TGFBI stands for Transforming Growth Factor Beta Induced.
It is a gene that gets switched on in response to **TGF-beta** signaling. My survival analysis showed that pancreatic
cancer patients with high TGFBI expression in their tumors survive for
shorter periods of time.

## Methodology

### Data Used
- TCGA-PAAD gene expression and survival data (183 patients) from UCSC Xena
- Annotated Seurat object from Person 2 (8344 cells, 5 cell types)

### Survival Analysis
I got TGFBI expression for all 183 patients and split them into  certain categories (high and
low). Then I drew a Kaplan-Meier curve to compare
survival between the two groups using R langauge.

### CellChat Analysis
I ran CellChat on the single cell data from Person 2 to find out which cells
are sending TGF-beta signals to Schwann cells.

### Trajectory Analysis
I took ran Monocle3 pseudotime analysis so basically  I ran clustering, learned the
trajectory graph via R, and ordered cells by pseudotime. This shows how Schwann cells change their state over time along a trajectory

## Results

### Survival Analysis
The Kaplan-Meier plot shows overall survival of 183 pancreatic cancer
patients split by TGFBI expression. Patients with high TGFBI (red line)
show slightly worse survival compared to low TGFBI (blue line). The
p-value is 0.083 which means the trend exists but is not strongly
significant — this could be because the sample size is small (183 patients).
This still supports the paper's finding that TGFBI is linked to worse outcomes.

### CellChat TGF-beta Signaling
The circle plot shows TGF-beta signaling between 5 cell types. Tumor cells
and Macrophages send the strongest TGF-beta signals (thick lines) to other
cells including Schwann cells. Schwann cells also signal back to themselves
(the loop). This confirms the paper's finding that TGF-beta signaling
from the tumor microenvironment reaches Schwann cells and likely drives
the TGFBI+ Schwann cell subtype.

### Monocle3 Trajectory
The UMAP plot shows Schwann cells colored by pseudotime (dark blue = early,
yellow = late). The trajectory line shows Schwann cells start from one
state (dark blue, top left) and gradually transition to a different state
(yellow, bottom). This suggests Schwann cells go through a developmental
change inside the tumor — possibly transitioning toward the TGFBI+ subtype
as pseudotime increases.
