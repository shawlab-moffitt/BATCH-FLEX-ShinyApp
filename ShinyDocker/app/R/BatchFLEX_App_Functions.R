# Function to remove genes that are lowly expressed
ExprFilter2 <- function(vec,criteria,proportion) {
  Samp2meet <- length(vec)*proportion
  meet <- sum(vec>criteria)
  if (meet > Samp2meet) {
    return(TRUE)
  } else { return(FALSE) }
}

cv <- function(x){
  (sd(x)/mean(x))*100
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

MouseToHuman <- function(data,conv,gene_col = TRUE) {
  require(dplyr)
  data <- as.data.frame(data)
  if (!gene_col) {
    data <- cbind(Genes = rownames(data),data)
  }
  data <- merge(conv,data,by.x = "Mouse", by.y = colnames(data)[1], all.y = T)
  data[which(is.na(data[,"Human"])),"Human"] <- toupper(data[which(is.na(data[,"Human"])),"Mouse"])
  data <- data[,-1]
  colnames(data)[1] <- "Genes"
  if (TRUE %in% duplicated(data[, 1])) {
    data_dup <- data %>% dplyr::group_by(Genes) %>% dplyr::filter(n() > 1) %>% as.data.frame()
    data_nodup <- data %>% dplyr::group_by(Genes) %>% dplyr::filter(n() == 1) %>% as.data.frame()
    data_dup_summ <- data_dup %>%
      dplyr::group_by(Genes) %>%
      dplyr::summarise_all(max) %>%
      as.data.frame()
    data <- rbind(data_nodup,data_dup_summ)
  }
  rownames(data) <- data[,1]
  data <- as.matrix(data)
  return(data)
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

detect_species <- function(genes) {
  # check the capitalization pattern
  is_human_gene <- function(gene) {
    return(grepl("^[A-Z]+$", gene))
  }
  is_mouse_gene <- function(gene) {
    return(grepl("^[A-Z][a-z]*$", gene))
  }
  # Count the number of human and mouse gene patterns
  human_count <- sum(sapply(genes, is_human_gene))
  mouse_count <- sum(sapply(genes, is_mouse_gene))
  # Determine the majority match
  if (human_count > mouse_count) {
    return("human")
  } else if (mouse_count > human_count) {
    return("mouse")
  } else {
    return("undetermined") # If counts are equal or if there are no matches
  }
}