#' buildGODatabase
#' 
#' @description driver to build GO database
#' 
#' @param goa character string path name to downloaded goa_human.gaf
#' @param gobasic character string path name to downloaded go-basic.obo
#' @param dir character string path name to directory to hold subdirectory GODB_RDATA
#' @param verbose Boolean if TRUE print out some diagnostic info
#' 
#' @details
#' download goa_human.gaf from https://current.geneontology.org/products/pages/downloads.html
#' download go-basic.obo from https://geneontology.org/docs/download-ontology/
#' parameter dir should be omitted or NULL except for the developer harvesting the updated .RData DBs
#' 
#' @examples
#' \dontrun{
#' # replace my path names for goa and gobasic with your own!!
#' # these were obtained from the download sites listed in 'details' section
#' goa<-"~/goa_human.gaf"
#' gobasic<-"~/go-basic.obo"
#' GOGOA<-buildGODatabase(goa,gobasic,dir="~/",verbose=FALSE)
#' # > dim(GOGOA)
#' # [1] 720139      5
#' # > GOGOA[1:5,]
#' #      HGNC          GO           RELATION      NAME                    ONTOLOGY            
#' # [1,] "NUDT4B"      "GO:0003723" "enables"     "RNA binding"           "molecular_function"
#' # [2,] "NUDT4B"      "GO:0005515" "enables"     "protein binding"       "molecular_function"
#' # [3,] "NUDT4B"      "GO:0046872" "enables"     "metal ion binding"     "molecular_function"
#' # [4,] "NUDT4B"      "GO:0005829" "located_in"  "cytosol"               "cellular_component"
#' # [5,] "TRBV20OR9-2" "GO:0002376" "involved_in" "immune system process" "biological_process"
#' }
#' 
#' # here is a small example that you can run
#' f1<-system.file("extdata","goa_human.small.gaf",package="minimalistGODB")
#' f2<-system.file("extdata","go-basic.small.obo",package="minimalistGODB")
#' GODBsmall<-buildGODatabase(f1,f2,verbose=FALSE)
#' 
#' @details
#' The output GOGOA was saved as an .RData file.
#' This was too large for CRAN.
#' It is available from https://github.com/barryzee/GO 
#' 
#' @return returns GO database with columns c("HGNC","GO","RELATION","NAME","ONTOLOGY")
#' 
#' @export
buildGODatabase<-
	function(goa,gobasic,dir=NULL,verbose=FALSE) {
		GOA<-parseGOA(goa)
		GO<-parseGOBASIC(gobasic,verbose)
		GOGOA<-joinGO(GOA,GO)
		GOGOA3<-subsetGOGOA(GOGOA)
		
		if(!is.null(dir)) {
		  subd<-sprintf("%s/%s",dir,"GODB_RDATA")
		  if(!dir.exists(subd))
		    dir.create(subd)
		  save(GOA,file=sprintf("%s/GOA.RData",subd))
		  save(GO,file=sprintf("%s/GO.RData",subd))
		  save(GOGOA,file=sprintf("%s/GOGOA.RData",subd))
		  save(GOGOA3,file=sprintf("%s/GOGOA3.RData",subd))
		}
		
		return(GOGOA3)
	}

#' parseGOA
#' 
#' @description parse goa_human.gaf
#' 
#' @param goa character string path name to downloaded goa_human.gaf
#' 
#' @details
#' download goa_human.gaf from https://current.geneontology.org/products/pages/downloads.html
#' 
#' @examples
#' \dontrun{
#' # replace my path name for goa with your own!!
#' # this was obtained from the download sites listed in 'details' section
#' GOA<-parseGOA("~/goa_human.gaf")
#' # GOA[1:5,]
#' #      HGNC          GO           RELATION     
#' # [1,] "NUDT4B"      "GO:0003723" "enables"    
#' # [2,] "NUDT4B"      "GO:0005515" "enables"    
#' # [3,] "NUDT4B"      "GO:0046872" "enables"    
#' # [4,] "NUDT4B"      "GO:0005829" "located_in" 
#' # [5,] "TRBV20OR9-2" "GO:0002376" "involved_in"
#' }
#' # here is a small example that you can run
#' f<-system.file("extdata","goa_human.small.gaf",package="minimalistGODB")
#' GOAsmall<-parseGOA(f)
#' 
#'
#' @return returns matrix with columns c("HGNC","GO","RELATION")
#' 
#' @export
parseGOA<-
	function(goa) {
		x<-strsplit(grep("^UniProtKB",readLines(goa),value=TRUE),"\t")
		lx<-length(x)
		m<-matrix(nrow=lx,ncol=3)
		colnames(m)<-c("HGNC","GO","RELATION")
		for(i in 1:lx) {
			m[i,"HGNC"]<-x[[i]][3]
			m[i,"GO"]<-x[[i]][5]
			m[i,"RELATION"]<-x[[i]][4]
		}
		return(unique(m))
	}

#' parseGOBASIC
#' 
#' @description parse go-basic.obo
#' 
#' @param gobasic character string path name to downloaded go-basic.obo
#' @param verbose Boolean if TRUE print out some diagnostic info
#' 
#' @details
#' download go-basic.obo from https://geneontology.org/docs/download-ontology/
#' 
#' @examples
#' \dontrun{
#' # replace my path name for gobasic with your own!!
#' # this was obtained from the download sites listed in 'details' section
#' GO<-parseGOBASIC("~/go-basic.obo",verbose=FALSE)
#' # GO$bp[1:5,]
#' #            GO           NAME                               ONTOLOGY            
#' # GO:0000001 "GO:0000001" "mitochondrion inheritance"        "biological_process"
#' # GO:0000002 "GO:0000002" "mitochondrial genome maintenance" "biological_process"
#' # GO:0000011 "GO:0000011" "vacuole inheritance"              "biological_process"
#' # GO:0000012 "GO:0000012" "single strand break repair"       "biological_process"
#' # GO:0000017 "GO:0000017" "alpha-glucoside transport"        "biological_process"
#' }
#'
#' # here is a small example that you can run
#' f<-system.file("extdata","go-basic.small.obo",package="minimalistGODB")
#' GOsmall<-parseGOBASIC(f)
#' 
#' @return returns a list whose components are c("m",  "bp", "mf", "cc")
#' 
#' @export
parseGOBASIC<-
	function(gobasic,verbose=FALSE) {
		l<-list()
		x<-readLines(gobasic)
		
		v<-grep("[Term]",x,fixed=TRUE)
		lv<-length(v)
		
		if(verbose)
			message(c("number of terms:",lv))
		
		m<-matrix(nrow=lv,ncol=3)
		colnames(m)<-c("GO","NAME","ONTOLOGY")
		for(i in 1:lv) {
			m[i,"GO"]<-x[v[i]+1]
			m[i,"NAME"]<-x[v[i]+2]
			m[i,"ONTOLOGY"]<-x[v[i]+3]
		}
		if(verbose) {
			message(c("intial m"))
			message(m[1:10,])
		}
		
		# filter out "obsolete" in "NAME"
		g<-grep("obsolete",m[,"NAME"])
		if(length(g)>0)
		  m<-m[-g,]
		if(verbose) {
			message(c("non obsolete m"))
			message(m[1:10,])
		}
		
		m[,"GO"]<-substring(m[,"GO"],5,14)
		m[,"NAME"]<-substring(m[,"NAME"],7,1000)
		m[,"ONTOLOGY"]<-substring(m[,"ONTOLOGY"],12,50)
		
		rownames(m)<-m[,"GO"]
			
		bp<-grep("biological_process",m[,"ONTOLOGY"])
		if(verbose) {
			message(c("BP indices",length(bp)))
			message(bp[1:10])
			message(m[bp[1:10],])
		}		
		
		l$m<-m
		l$bp<-m[bp,]
		mf<-grep("molecular_function",m[,"ONTOLOGY"])
		l$mf<-m[mf,]
		cc<-grep("cellular_component",m[,"ONTOLOGY"])
		l$cc<-m[cc,]	
		
		# https://geneontology.org/stats.html
		if(verbose) {
			message("THEORETICAL NUMBER OF TERMS 26091 10154 4022")
			lx<-length(x)
			message(c("LX",lx,"LBP",length(bp),"LMF",length(mf),"LCC",length(cc)))
		}
		
		return(l)
	}

#' joinGO
#' 
#' @description join the outputs of parseGOA and parseGOBASIC to add the GO category name and the ontology to GOA
#' 
#' @param GOA output of parseGOA()
#' @param GO output of parseGOBASIC() 
#' 
#' @examples
#' GOGOA<-joinGO(GOA,GO)
#' # GOGOA[1:5,]
#' # HGNC          GO           RELATION      NAME                    ONTOLOGY            
#' # [1,] "NUDT4B"      "GO:0003723" "enables"     "RNA binding"           "molecular_function"
#' # [2,] "NUDT4B"      "GO:0005515" "enables"     "protein binding"       "molecular_function"
#' # [3,] "NUDT4B"      "GO:0046872" "enables"     "metal ion binding"     "molecular_function"
#' # [4,] "NUDT4B"      "GO:0005829" "located_in"  "cytosol"               "cellular_component"
#' # v[5,] "TRBV20OR9-2" "GO:0002376" "involved_in" "immune system process" "biological_process"
#' # GO_NAME                              
#' # [1,] "GO_0003723__RNA_binding"          
#' # [2,] "GO_0005515__protein_binding"      
#' # [3,] "GO_0046872__metal_ion_binding"    
#' # [4,] "GO_0005829__cytosol"              
#' # [5,] "GO_0002376__immune_system_process"
#' 
#' # querying GOGOA to compute gene enrichment of some GO categories
#' hgncList<-GOGOA[1:1000,"HGNC"]
#' ontology<-"biological_process"
#' w<-which(GOGOA[,"ONTOLOGY"] == ontology)
#' GOGOA<-GOGOA[w,]
#' w<-which(GOGOA[,"HGNC"] %in% hgncList)
#' t<-sort(table(GOGOA[w,"NAME"]),decreasing=TRUE)[1:10]
#'
#' @return returns a matrix with columns c("HGNC","GO","RELATION","NAME","ONTOLOGY")
#' 
#' @export
joinGO<-
	function(GOA,GO) {
		m<-matrix(nrow=nrow(GOA),ncol=5)
		NAME<-GO$m[GOA[,"GO"],"NAME"]
		gx<-cbind(GOA,NAME)
		rownames(gx)<-NULL
		ONTOLOGY<-GO$m[GOA[,"GO"],"ONTOLOGY"]
		gxx<-cbind(gx,ONTOLOGY)
		rownames(gxx)<-NULL
		
		# generate 'safe' (suitable to be used as a variable or file name) combined GO category name
		#GO_NAME<-sprintf("%s__%s",gsub(":","_",gxx[,"GO"]),gsub(" ","_",gxx[,"NAME"]))
		GO_NAME<-sprintf("%s__%s",gsub(":","_",gxx[,"GO"]),gsub("[ /,]","_",gxx[,"NAME"]))
		gxxx<-cbind(gxx,GO_NAME)
		
		return(gxxx)
	}

#' subsetGOGOA
#' 
#' @description split GOGOA into 3 separate ontologies
#' 
#' @param GOGOA return value of minimalistGODB::joinGO()
#' 
#' @examples
#' #load("data/GOGOAsmall.RData")
#' GOGOA3small<-subsetGOGOA(GOGOAsmall)
#' 
#' @return returns a list containing subsets of GOGOA for each ontology, unique gene and cat lists, and stats
#' 
#' @export
subsetGOGOA<-
  function(GOGOA) {
    l<-list()
    
    l$ontologies<-list()
    ontologies<-unique(GOGOA[,"ONTOLOGY"])
    for(ONTOLOGY in ontologies) {
      w<-which(GOGOA[,"ONTOLOGY"] == ONTOLOGY)
      l$ontologies[[ONTOLOGY]]<-as.matrix(GOGOA[w,c("HGNC","GO","RELATION","NAME","GO_NAME")])
    }
    
    l$cats<-list()
    l$genes<-list()
    l$stats<-list()
    l$stats$ncats<-list()
    l$stats$ngenes<-list()
    for(ONTOLOGY in ontologies) {
      l$cats[[ONTOLOGY]]<-unique(l$ontologies[[ONTOLOGY]][,"GO_NAME"])
      l$stats$ncats[[ONTOLOGY]]<-length(l$cats[[ONTOLOGY]])
      l$genes[[ONTOLOGY]]<-unique(l$ontologies[[ONTOLOGY]][,"HGNC"])
      l$stats$ngenes[[ONTOLOGY]]<-length(l$genes[[ONTOLOGY]])
    }
    l$genes[["ALL"]]<-unique(GOGOA[,"HGNC"])
    l$cats[["ALL"]]<-unique(GOGOA[,"GO_NAME"])
    return(l)
  }