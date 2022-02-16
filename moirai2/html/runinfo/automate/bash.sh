# Retrieve runinfo and summary from accessionIDs
mkdir -p runinfo/output

perl moirai2.pl \
-c edirect.sif \
-i '$accessionid->studyid->$studyid' \
-o '$accessionid->runinfo->$runinfo,$accessionid->summaryfile->$summary' \
command/edirect/studyid_to_sra_info_files.json \
'$runinfo=runinfo/output/$accessionid.runinfo.txt' \
'$summary=runinfo/output/$accessionid.summary.txt'
