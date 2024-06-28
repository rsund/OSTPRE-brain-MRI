library(dplyr)

basepath <- "/research/groups/ostpre/ostpre_bids/"

tree <- readr::read_csv(paste0(basepath, "dicom/tree"), col_names=FALSE) |>
  mutate(
    num=stringr::str_remove_all(X1,"^\\./sub-|/ses|:"),
  ) |>
  tidyr::separate(num, c("lomno1", "ses"), "-") |>
  mutate(
    lomno1=as.numeric(lomno1),
    ses=as.numeric(ses),
    inpath=paste0(basepath, "dicom/sub-", lomno1, "/ses-", ses),
    subpath=paste0(basepath, "nifti/sub-", sprintf('%05d', lomno1)),
    sespath=paste0(subpath, "/ses-", sprintf('%02d', ses))
  ) |>
  arrange(lomno1, ses)

tags <- NULL

for (i in 1:nrow(tree)) { 

  lid <- as.numeric(tree[i, "lomno1"])
  sid <- as.numeric(tree[i, "ses"])
  supa <- as.character(tree[i, "subpath"])
  if (!dir.exists(supa)) dir.create(supa)
    
  inpa <- as.character(tree[i,"inpath"])
  oupa <- as.character(tree[i,"sespath"])
  cmd <- paste("~/dcm2niix -z y -ba n -o", oupa, inpa, ">>", paste0(oupa,"/log.txt"), "2>&1")
  
  unlink(oupa, recursive = TRUE)
  dir.create(oupa)
  system(cmd)
  
  loki <- grep("^Convert",readLines(paste0(oupa, "/log.txt")),value=TRUE)
  jsonfiles <- Sys.glob(paste0(oupa, "/*.json"))
  
  for(jsonfil in jsonfiles) {
    
    niifil <- gsub(".json$",".nii.gz",jsonfil)
    filesize <- file.size(niifil)
    dime <- stringr::word(grep(gsub(".json$"," ",jsonfil), loki, value=TRUE, fixed=TRUE), -1)
    
    tagstmp <- as.list(jsonlite::fromJSON(txt=jsonfil)) |>
      purrr::map_if(is.vector, list) |>
      tibble::as_tibble() |>
      mutate(
        Path=gsub(".*nifti/", "", jsonfil),
        lomno1=lid,
        ses=sid,
        NiftiPath=gsub(".*nifti/", "", niifil),
        FileSize=filesize,
        ConvertDimensions=dime
        )
    
    tags <- tags |>
      union_all(tagstmp)
  }
  
}

qs::qsavem(tags, file=paste0(basepath, "nifti/tags.qs"))
openxlsx::write.xlsx(tags,"../nifti/tags.xlsx")
