# wheel-running-analysis
MATLAB scripts for analyzing TDT fiber-photometry data.

"TDTbin2mat.m", "TDTfilter.m", and "SEV2mat.m" must be in the same directory as the analysis scripts to work.

Instructions for each analysis file are included as comments within the scripts.

BATCH PROCESSING:
Batch processing allows you to analyze multiple tanks within a folder. To keep the data output organized, it is best to batch process folders
with smaller batches of tanks. For example, batch process all the "Locked" tanks in the "Post_Coc" parent folder followed by the "Unlocked." Repeat
for each parent tank folder.

To use batch processing, the recommended folder organization is as follows:
![Folder Organization for Batch Processing](/img/Batch_Folder_Org.png)

