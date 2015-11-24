### R code from vignette source 'NewSchema.Rnw'

###################################################
### code chunk number 1: NewSchema.Rnw:27-29
###################################################
options(width=60)
options(continue="  ")


###################################################
### code chunk number 2: load packages
###################################################
library(AnnotationDbi)


###################################################
### code chunk number 3: create dir
###################################################
package <- "targetscan.Hs.eg.db"
packagedir <- paste(tempdir(), sep=.Platform$file.sep, "targetscan")
unlink(packagedir, recursive=TRUE)


###################################################
### code chunk number 4: copy template
###################################################
start.template <- system.file("AnnDbPkg-templates/KEGG.DB", 
                              package="AnnotationDbi")
tmp.keggdir <- paste(tempdir(), sep=.Platform$file.sep, "KEGG.DB")
unlink(tmp.keggdir, recursive=TRUE)
file.copy(start.template, tempdir(), recursive=TRUE)
file.rename(tmp.keggdir, packagedir)


###################################################
### code chunk number 5: update DESCRIPTION file
###################################################
desc <- packageDescription("targetscan", tempdir())
desc$License <- "BSD"
desc$Description <- 
  "@PKGTITLE@ assembled using data from TargetScan"
write.dcf(unclass(desc), attr(desc, "file"))


###################################################
### code chunk number 6: bimap code
###################################################
zzzfile <- paste(packagedir, sep=.Platform$file.sep, 
                 "R", "zzz.R")
bimap.code <- '
### Mandatory fields: objName, Class and L2Rchain
TARGETSCAN_DB_AnnDbBimap_seeds <- list(
    list(objName="MIRBASE2FAMILY",
         Class="AnnDbBimap",
         L2Rchain=list(
           list(tablename="mirbase",
                Lcolname="mirbase_id",
                Rcolname="family"
                ),
           list(
                tablename="mirna_family",
                Lcolname="_id",
                Rcolname="name"
                )
           )
         ),
     list(objName="MIRNA",
          Class="miRNAAnnDbBimap",
          L2Rchain=list(
            list(tablename="mirbase",
                 Lcolname="mirbase_id",
                 Rcolname="mirbase_id",
                 Rattribnames=c(
                   MiRBase_Accession="{mirbase_accession}",
                   Seed_m8="{seed_m8}",
                   Species_ID="species.name",
                   Mature_sequence="{mature_sequence}",
                   Family_conservation="{family_conservation}"),
                 Rattrib_join=
                   "LEFT JOIN species ON {species}=species.id"
                 )
            )
          ),
      list(objName="TARGETS",
           Class="AnnDbBimap",
           L2Rchain=list(
             list(tablename="genes",
                  Lcolname="gene_id",
                  Rcolname="_id"
                  ),
             list(tablename="targets",
                  Lcolname="target",
                  Rcolname="family"
                  ),
             list(tablename="mirna_family",
                  Lcolname="_id",
                  Rcolname="name"
                  )
             )
           ),
      list(objName="TARGETSFULL",
           Class="miRNATargetAnnDbBimap",
           L2Rchain=list(
             list(tablename="genes",
                  Lcolname="gene_id",
                  Rcolname="_id"
                  ),
             list(tablename="targets",
                  Lcolname="target",
                  Rattribnames=c(
                    UTR_start="{utr_start}",
                    UTR_end="{utr_end}",
                    MSA_start="{msa_start}",
                    MSA_end="{msa_end}",
                    Seed_match="seed_match.name",
                    PCT="{pct}"),
                  Rattrib_join=paste(sep="", 
                    "LEFT JOIN seed_match ON {seed_match}=seed_match._id ",
                    "LEFT JOIN mirna_family AS _R ON {family}=_R._id"),
                  Rcolname="name"
##                   ),
##              list(tablename="mirna_family",
##                   Lcolname="_id",
##                   Rcolname="name"
                  )
             )
           )
)
'
cat(bimap.code, file=zzzfile, append=TRUE)


###################################################
### code chunk number 7: add code to create ann objects at load time
###################################################
create.code <- '
createAnnObjs.TARGETSCAN_DB <- function(prefix, objTarget, 
                                        dbconn, datacache)
{
    ## checkDBSCHEMA(dbconn, "TARGETSCAN_DB")

    ## AnnDbBimap objects
    seed0 <- list(
        objTarget=objTarget,
        datacache=datacache
    )
    ann_objs <- AnnotationDbi:::createAnnDbBimaps(
                  TARGETSCAN_DB_AnnDbBimap_seeds, seed0)

    ## Reverse maps
    revmap2 <- function(from, to)
    {
        map <- revmap(ann_objs[[from]], objName=to)
        L2Rchain <- map@L2Rchain
        tmp <- L2Rchain[[1]]@filter
        L2Rchain[[1]]@filter <-
          L2Rchain[[length(L2Rchain)]]@filter
        L2Rchain[[length(L2Rchain)]]@filter <- tmp
        map@L2Rchain <- L2Rchain
        map
    }
    ann_objs$FAMILY2MIRBASE <- revmap2("MIRBASE2FAMILY", 
                                       "FAMILY2MIRBASE")
    
    ## 1 special map that is not an AnnDbBimap object 
    ## (just a named integer vector)
    ann_objs$MAPCOUNTS <- 
      AnnotationDbi:::createMAPCOUNTS(dbconn, prefix)

    AnnotationDbi:::prefixAnnObjNames(ann_objs, prefix)
}
'
cat(create.code, file=zzzfile, append=TRUE)


###################################################
### code chunk number 8: remove createAnnObjs.SchmemaChoice function
###################################################
zzz <- readLines(zzzfile)
zzz <- sub('createAnnObjs.SchemaChoice("@DBSCHEMA@",', 
           'createAnnObjs.@DBSCHEMA@(',
           zzz, fixed=TRUE)
writeLines(zzz, zzzfile)


###################################################
### code chunk number 9: implement new bimap classes
###################################################
class.code <- '
setClass("miRNAAnnDbBimap", contains="AnnDbBimap")
setClass("miRNATargetAnnDbBimap", contains="AnnDbBimap")

setMethod("as.list", "miRNATargetAnnDbBimap",
          function(x, ...)
          {
            y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
            makemiRNATargetNode <- function(name, UTR_start, UTR_end,
                                            MSA_start, MSA_end, 
                                            Seed_match, PCT, ...)
              {
                l <- mapply(list, SIMPLIFY=FALSE,
                            miR.Family=name,
                            UTR.start=UTR_start,
                            UTR.end=UTR_end,
                            MSA.start=MSA_start,
                            MSA.end=MSA_end,
                            Seed.match=Seed_match,
                            PCT=PCT)
                for (i in seq_along(l)) { 
                  class(l[[i]]) <- "targetscan.Target" 
                }
                l
              }
            AnnotationDbi:::.toListOfLists(y, mode=1, 
                                           makemiRNATargetNode)
          }
          )

setMethod("as.list", "miRNAAnnDbBimap",
          function(x, ...)
          {
            y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
            makemiRNANode <- function(mirbase_id, 
                                      MiRBase_Accession,
                                      Seed_m8, Species_ID, 
                                      Mature_sequence,
                                      Family_conservation, ...)
              {
                l <- list(MiRBase.ID=mirbase_id,
                          MiRBase.Accession=MiRBase_Accession,
                          Seed.m8=Seed_m8,
                          Species=Species_ID,
                          Mature.sequence=Mature_sequence,
                          Family.conservation=Family_conservation,
                          ...)
                class(l) <- "targetscan.MiRBase"
                l
              }
            AnnotationDbi:::.toListOfLists(y, mode=1, makemiRNANode)
          }
          )
'
zzz <- readLines(zzzfile)
insert.class <- grep('require[(]"methods', zzz)[1]
zzz <- c(zzz[1:insert.class], 
         paste("   ", strsplit(class.code, "\n")[[1]]),
         zzz[(insert.class+1):length(zzz)])
writeLines(zzz, zzzfile)


###################################################
### code chunk number 10: add print methods
###################################################
print.code <- '
print.targetscan.MiRBase <- function(x, ...) {
  for (name in names(x)) {
    name2 <- sub("\\\\.", " ", name)
    s <- paste(sep="", name2, ": ", x[[name]])
    cat(strwrap(s, exdent=4), sep="\n")
  }
  invisible(x)
}

print.targetscan.Target <- function(x, ...) {
  for (name in names(x)) {
    name2 <- sub("\\\\.", " ", name)
    s <- paste(sep="", name2, ": ", x[[name]])
    cat(strwrap(s, exdent=4), sep="\n")
  }
  invisible(x)
}
'
cat(print.code, file=zzzfile, append=TRUE)


###################################################
### code chunk number 11: update NAMESPACE file
###################################################
namespace.file <- zzzfile <- paste(packagedir, sep=.Platform$file.sep, 
                                   "NAMESPACE")
namespace.text <- '
exportClasses("miRNAAnnDbBimap",
              "miRNATargetAnnDbBimap")

export(print.targetscan.MiRBase, print.targetscan.Target)
S3method("print", targetscan.MiRBase)
S3method("print", targetscan.Target)
'
cat(namespace.text, file=namespace.file, append=TRUE)


###################################################
### code chunk number 12: remove manual pages
###################################################
mandir <- paste(tempdir(), sep=.Platform$file.sep, 
                "targetscan", "man")
unlink(mandir, recursive=TRUE)


###################################################
### code chunk number 13: add database schema
###################################################
schema.text <- '
--
-- TARGETSCAN_DB schema
-- ====================
--

CREATE TABLE genes (
  _id INTEGER PRIMARY KEY,
  gene_id VARCHAR(10) NOT NULL UNIQUE          -- Entrez Gene ID
);

CREATE TABLE mirna_family (
  _id INTEGER PRIMARY KEY,
  name VARCHAR(255) NOT NULL                   -- miRNA family
);

CREATE TABLE species (
  id VARCHAR(10) PRIMARY KEY,                  -- species ID from NCBI
  name VARCHAR(100) NOT NULL                   -- species name
);

CREATE TABLE seed_match (
  _id INTEGER PRIMARY KEY,
  name VARCHAR(10)
);

CREATE TABLE mirbase (
  mirbase_id VARCHAR(50) PRIMARY KEY,          -- MiRBase ID
  mirbase_accession CHAR(12) NOT NULL UNIQUE,  -- MiRBase accession
  family INTEGER NOT NULL,                     -- REFERENCES family
  seed_m8 CHAR(7) NOT NULL,                    -- seed m8
  species VARCHAR(10) NOT NULL,                -- REFERENCES species
  mature_sequence VARCHAR(100) NOT NULL,       -- mature sequence
  family_conservation VARCHAR(3) NOT NULL,     -- family convervation
  FOREIGN KEY (family) REFERENCES mirna_family (_id),
  FOREIGN KEY (species) REFERENCES species (id)
);

CREATE TABLE targets (
  family INTEGER NOT NULL,                     -- REFERENCES family
  target INTEGER NOT NULL,                     -- REFERENCES genes
  species VARCHAR(10) NOT NULL,                -- REFERENCES species
  utr_start INTEGER NOT NULL,                  -- UTR start
  utr_end INTEGER NOT NULL,                    -- UTR end
  msa_start iNTEGER NOT NULL,                  -- MSA start
  msa_end INTEGER NOT NULL,                    -- MSA end
  seed_match INTEGER NOT NULL,                 -- REFERENCES seed_match
  pct VARCHAR(10) NOT NULL,                    -- PCT
  FOREIGN KEY (family) REFERENCES mirna_family (_id),
  FOREIGN KEY (target) REFERENCES genes (_id),
  FOREIGN KEY (species) REFERENCES species (id),
  FOREIGN KEY (seed_match) REFERENCES seed_match (_id)
);

-- Metadata tables.
CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);
CREATE TABLE map_counts (
  map_name VARCHAR(80) PRIMARY KEY,
  count INTEGER NOT NULL
);
CREATE TABLE map_metadata (
  map_name VARCHAR(80) NOT NULL,
  source_name VARCHAR(80) NOT NULL,
  source_url VARCHAR(255) NOT NULL,
  source_date VARCHAR(20) NOT NULL
);
-- Indexes
'
dir.create(paste(packagedir, sep=.Platform$file.sep, 
                 "inst", "DBschemas", "schemas_1.0"),
           recursive=TRUE, mode="0755")
cat(schema.text, file=paste(packagedir, sep=.Platform$file.sep,
                   "inst", "DBschemas", "schemas_1.0", 
                   "TARGETSCAN_DB.sql"))


###################################################
### code chunk number 14: place of database file
###################################################
dbfile <- file.path(tempdir(), "targetscan.Hs.eg.sqlite")
unlink(dbfile)


###################################################
### code chunk number 15: download data (eval = FALSE)
###################################################
## ## Download TargetScan data
## targetfile <- file.path(tempdir(), "Predicted_Targets_Info.txt")
## targetfile.zip <- paste(targetfile, sep="", ".zip")
## if (!file.exists(targetfile)) {
##   data.url <- paste(sep="", "http://www.targetscan.org/vert_50/",
##          "vert_50_data_download/Predicted_Targets_Info.txt.zip")
##   download.file(data.url, destfile=targetfile.zip)
##   targetfile.tmp <- 
##     zip.file.extract(targetfile, basename(targetfile.zip))
##   file.copy(targetfile.tmp, targetfile)
## }
## 
## familyfile <- file.path(tempdir(), "miR_Family_Info.txt")
## familyfile.zip <- paste(familyfile, sep="", ".zip")
## if (!file.exists(familyfile)) {
##   data.url <- paste(sep="", "http://www.targetscan.org/vert_50/",
##                     "vert_50_data_download/miR_Family_Info.txt.zip")
##   download.file(data.url, destfile=familyfile.zip)
##   familyfile.tmp <- 
##     zip.file.extract(familyfile, basename(familyfile.zip))
##   file.copy(familyfile.tmp, familyfile)
## }
## 
## taxfile <- file.path(tempdir(), "names.dmp")
## taxfile.zip <- file.path(tempdir(), "taxdmp.zip")
## if (!file.exists(taxfile)) {
##   data.url <- "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
##   download.file(data.url, destfile=taxfile.zip)
##   taxfile.tmp <- 
##     zip.file.extract(taxfile, basename(taxfile.zip))
##   file.copy(taxfile.tmp, taxfile)
## }


###################################################
### code chunk number 16: read in data (eval = FALSE)
###################################################
## family <- read.delim(familyfile)
## targets <- read.delim(targetfile)
## 
## tax <- read.delim(taxfile, header=FALSE, quote="")
## tax <- tax[,-c(2,4,6,8)]
## species <- unique(family$Species.ID)
## names <- tax[,2][ match(species, tax[,1]) ]
## species <- data.frame(id=species, name=names)


###################################################
### code chunk number 17: create database and tables
###################################################
## Create the database file
library(RSQLite)
drv <- dbDriver("SQLite")
db <- dbConnect(drv, dbname=dbfile)

## Create tables
create.sql <- strsplit(schema.text, "\n")[[1]]
create.sql <- paste(collapse="\n", create.sql)
create.sql <- strsplit(create.sql, ";")[[1]]
create.sql <- create.sql[-length(create.sql)] # nothing to run here
tmp <- sapply(create.sql, function(x) sqliteQuickSQL(db, x))


###################################################
### code chunk number 18: put the data into the tables (eval = FALSE)
###################################################
## ## Append data
## ## Species
## dbBeginTransaction(db)
## dbGetPreparedQuery(db, 'INSERT INTO "species" VALUES(:id, :name);',
##                    bind.data=species)
## dbCommit(db)
## 
## ## miRNA families
## family$miR.family <- sub("^miR-141/200$", "miR-141/200a", 
##                          family$miR.family, fixed=FALSE)
## family$miR.family <- sub("^miR-200bc/420$", "miR-200bc/429",
##                          family$miR.family, fixed=FALSE)
## fam <- unique(as.character(family$miR.family))
## fam <- cbind(id=seq_along(fam), fam)
## dbBeginTransaction(db)
## dbGetPreparedQuery(db, 
##                    "INSERT INTO 'mirna_family' VALUES(:id, :fam);",
##                    bind.data=as.data.frame(fam))
## dbCommit(db)
## 
## ## mirbase table
## mirbase <- family[,c("MiRBase.ID", "MiRBase.Accession",
##                      "miR.family", "Seed.m8", "Species.ID", 
##                      "Mature.sequence",
##                      "Family.Conservation.")]
## mirbase$miR.family <- fam[,1][ match(family$miR.family, fam[,2]) ]
## colnames(mirbase) <- letters[1:7]
## dbBeginTransaction(db)
## dbGetPreparedQuery(db, bind.data=mirbase,
##                    'INSERT INTO "mirbase" VALUES(:a,:b,:c,:d,:e,:f,:g)')
## dbCommit(db)
## 
## ## keep entries for human only
## targets2 <- targets
## targets2 <- targets2[ targets2$Species.ID == 9606, ]
## targets2$Gene.ID[ targets2$Gene.ID == 934] <- 100133941
## 
## ## genes
## gs <- unique(targets2$Gene.ID)
## gs <- cbind(seq_along(gs), gs)
## dbBeginTransaction(db)
## dbGetPreparedQuery(db, bind.data=data.frame(a=gs[,1],
##                          b=as.integer(gs[,2])),
##                    "INSERT INTO 'genes' VALUES(:a,:b);")
## dbCommit(db)
## 
## ## seed_match
## sm <- sort(unique(as.character(targets$Seed.match)))
## sm <- cbind(seq_along(sm), sm)
## dbBeginTransaction(db)
## dbGetPreparedQuery(db, bind.data=data.frame(a=sm[,1], b=sm[,2]),
##                    "INSERT INTO 'seed_match' VALUES(:a,:b);")
## dbCommit(db)
## 
## ## targets, human only :(
## targets2$miR.Family <- fam[,1][ match(targets2$miR.Family, fam[,2]) ]
## targets2$Gene.ID <- gs[,1][ match(targets2$Gene.ID, gs[,2]) ]
## targets2$Seed.match <- sm[,1][ match(targets2$Seed.match, sm[,2]) ]
## colnames(targets2) <- sub(".", "_", colnames(targets2), fixed=TRUE)
## dbBeginTransaction(db)
## dbGetPreparedQuery(db, bind.data=targets2,
##                    paste(sep="", 
##                          "INSERT INTO targets VALUES(:miR_Family,",
##                          ":Gene_ID, :Species_ID,",
##                          ":UTR_start, :UTR_end,",
##                          ":MSA_start, :MSA_end,",
##                          ":Seed_match, :PCT);"))
## dbCommit(db)


###################################################
### code chunk number 19: append the metadata (eval = FALSE)
###################################################
## ## metadata
## metadata <- rbind(c("DBSCHEMA", "TARGETSCAN_DB"),
##                   c("ORGANISM", "Homo sapiens"),
##                   c("SPECIES", "Human"),
##                   c("DBSCHEMAVERSION", "1.0"),
##                   c("TARGETSCANSOURCENAME", "TargetScan"),
##                   c("TARGETSCANSOURCEURL",
##                     paste(sep="", "http://www.targetscan.org/cgi-bin/",
##                           "targetscan/data_download.cgi?db=vert_50")),
##                   c("TARGETSCANSOURCEDATE", 
##                     format(Sys.time(), "%Y-%b%d")),
##                   c("TARGETSCANVERSION", "5.0"))
## q <- paste(sep="", "INSERT INTO 'metadata' VALUES('", metadata[,1],
##            "','", metadata[,2], "');")
## tmp <- sapply(q, function(x) sqliteQuickSQL(db, x))
## 
## ## map_counts
## map.counts <- rbind(c("FAMILY2MIRBASE", nrow(fam)),
##                     c("MIRBASE2FAMILY", nrow(mirbase)),
##                     c("MIRNA", nrow(mirbase)),
##                     c("TARGETS", nrow(gs)),
##                     c("TARGETSFULL", nrow(gs))
##                     )
## q <- paste(sep="", "INSERT INTO 'map_counts' VALUES('", map.counts[,1],
##            "',", map.counts[,2], ");")
## tmp <- sapply(q, function(x) sqliteQuickSQL(db, x))


###################################################
### code chunk number 20: quick checks (eval = FALSE)
###################################################
## if (dbGetQuery(db, "SELECT COUNT(*) FROM species") != nrow(species)) {
##   stop("FOOBAR")
## }
## if (dbGetQuery(db, "SELECT COUNT(*) FROM mirna_family") != nrow(fam)) {
##   stop("FOOBAR")
## }
## if (dbGetQuery(db, "SELECT COUNT(*) FROM mirbase") != nrow(mirbase)) {
##   stop("FOOBAR")
## }
## if (dbGetQuery(db, "SELECT COUNT(*) FROM genes") != nrow(gs)) {
##   stop("FOOBAR")
## }
## if (dbGetQuery(db, "SELECT COUNT(*) FROM seed_match") != nrow(sm)) {
##   stop("FOOBAR")
## }
## if (dbGetQuery(db, "SELECT COUNT(*) FROM targets") != nrow(targets2)) {
##   stop("FOOBAR")
## }


###################################################
### code chunk number 21: disconnect
###################################################
## Disconnect
dbGetQuery(db, "VACUUM")
dbDisconnect(db)


###################################################
### code chunk number 22: build new package (eval = FALSE)
###################################################
## seed <- new("AnnDbPkgSeed",
##             Package = package,
##             Version = "5.0-1",
##             PkgTemplate = packagedir,
##             AnnObjPrefix = "targetscan.Hs.eg",
##             Title = "TargetScan miRNA target predictions for human",
##             Author = "Gabor Csardi <Gabor.Csardi@foo.bar>",
##             Maintainer = "Gabor Csardi <Gabor.Csardi@foo.bar>",
##             organism = "Homo sapiens",
##             species = "Human",
##             biocViews = "AnnotationData, FunctionalAnnotation",
##             DBschema = "TARGETSCAN_DB",
##             AnnObjTarget = "TargetScan (Human)",
##             manufacturer = "None",
##             manufacturerUrl = "None"
##             )
## 
## unlink(paste(tempdir(), sep=.Platform$file.sep, package), recursive=TRUE)
## makeAnnDbPkg(seed, dbfile, dest_dir = tempdir())


###################################################
### code chunk number 23: install new package (eval = FALSE)
###################################################
## install.packages(paste(tempdir(), sep=.Platform$file.sep, "targetscan.Hs.eg.db"), 
##                  repos=NULL)


###################################################
### code chunk number 24: try new package (eval = FALSE)
###################################################
## library(targetscan.Hs.eg.db)
## 
## mget(c("346389", "54715"), targetscan.Hs.egTARGETS)
## mget("miR-328", revmap(targetscan.Hs.egTARGETS))
## mget("346389", targetscan.Hs.egTARGETSFULL)
## toTable(targetscan.Hs.egTARGETS)[1:5,]
## toTable(targetscan.Hs.egTARGETSFULL)[1:5,]
## 
## mget("hsa-let-7a", targetscan.Hs.egMIRNA)
## 
## checkMAPCOUNTS("targetscan.Hs.eg.db")


###################################################
### code chunk number 25: SessionInfo
###################################################
toLatex(sessionInfo())


