### FDB_beta code in 3 dimensions

- This code is capable of finding the minimum electric field for a target barrier (in kcal/mol)
- The author recommends reading first the information about the input before jumping into the main script
- If the user wishes to swap between the work modes (With and without electron transfer processes), define the process as either reduction or oxidation and set the potential to 0. On top of that, at line 40 either delete the potential from the definition of E_fdb or move it to the next line with an intro. Finally, comment the lines 130 and 131
- The warnings due to the COMMONs can be ignored
