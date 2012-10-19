require("RCurl")
require("XML")


#search.keywords 
query.term=c("(wetland OR fen OR marsh) AND 16S rRNA AND (soil OR sediment) AND Japan")

# Search and fetch XML from GenBank
searchgenbank <- function(query.term) {
  # change spaces to + in query
  query.gsub <- gsub(" ", "+", query.term)
  # change single-quotes to URL-friendly %22
  query.gsub <- gsub("'","%22", query.gsub)
  # Perform search and save history, this will save PMIDS in history
  gen.esearch <- getURL(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=Nucleotide&term=",
                              query.gsub, "&usehistory=y", sep = ""))
  # Parse esearch XML
  gen.esearch <- xmlTreeParse(gen.esearch, asText = TRUE)
  # Count number of hits (super assign)
  gen.count <- as.numeric(xmlValue(gen.esearch[["doc"]][["eSearchResult"]][["Count"]]))
  # Save WebEnv-string, it contains "links" to all articles in my search
  webenv<- xmlValue(gen.esearch[["doc"]][["eSearchResult"]][["WebEnv"]])
  key<- xmlValue(gen.esearch[["doc"]][["eSearchResult"]][["QueryKey"]])
  # Show how many articles that's being downloaded
  cat("Searching (downloading", gen.count, "articles)/n")
  
  ## We need to batch download, since efetch will cap at 10k articles ##
  # Start at 0
  RetStart<-0
  # End at 10k
  RetMax <-500
  # Calculate how many itterations will be needed
  Runs<- (gen.count %/%500) + 1
  # Create empty object
  gen.efetch <- NULL
  pub=NULL
  # Download XML based on hits saved in gen.esearch (WebEnv)
  for (i in 1:Runs) {
    # Save WebEnv-string, it contains "links" to all articles in my search
    x <- getURL(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&query_key=1&WebEnv=",webenv,"&retmode=xml&retstart=",RetStart,"&retmax=",500, sep =""))
    xml.data <- xmlTreeParse(x, useInternalNodes = TRUE)
    # Use xpathSApply to extract Journal name
    sequence=xpathSApply(xml.data,"//GBSet/GBSeq/GBSeq_sequence", xmlValue)
    #gb=xpathSApply(xml.data,"//GBSet/GBSeq/GBSeq_other-seqids/GBSeqid[1]", xmlValue)
    #gi=xpathSApply(xml.data,"//GBSet/GBSeq/GBSeq_other-seqids/GBSeqid[2]", xmlValue)
    gen.efetch <- c(gen.efetch, sequence)
    # Increase range for next batch
    RetStart<-RetStart+500
    
    cat("Now is the", i, "cycles./n")
  }
  # Print that download is completed
  cat("Completed download from genbank./n")
  
  # Return XML
  
  return(gen.efetch)
}


gen.esearch <- getURL(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=Nucleotide&term=",
                            query.gsub, "&usehistory=y", sep =""))

wetland=searchgenbank(query.term)
