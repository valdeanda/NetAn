
install.packages("visNetwork")
#devtools::install_github("datastorm-open/visNetwork") #for development version
library(visNetwork)

nodes<-read.table("../data/nodes.txt", sep="\t",header=TRUE)
edges<-read.table("../data/edges.txt", header = TRUE,sep="\t")



visNetwork(nodes,edges, height = "700px", width = "100%",main = "Test", submain = "Only 0.05% of the total interactions"
) %>%
  visEdges(arrows="to")%>%
  visOptions(selectedBy = "Domain", 
             highlightNearest = TRUE, 
             nodesIdSelection = TRUE,
             manipulation=TRUE) %>%
  visPhysics(stabilization =FALSE)%>%
visConfigure(enabled=TRUE)
