ped<-read.table("HLA-alleles-in-Ichip-cases.ped", sep="\t")
map<-read.table("HLA-alleles-in-Ichip-cases.map")

buildMapWithNames <- function(pedfile, mapfile){
        ped <- read.table(pedfile, sep = '\t', stringsAsFactors = FALSE, na.strings = c('0 0', ' ', NA, 'NA'))
        map <- read.table(mapfile, sep = '\t', stringsAsFactors = FALSE)
        names(ped)[7:ncol(ped)] <- as.character(map[, 2])
        return(ped)
}



convertAlleleStatus <- function(x){
        if(x %in% c('P P', 'A P', 'A A')){
                switch(x, 'P P' = 0, 'A P' = 1, 'A A' = 2)
        }
        else{
                return(x)
        }
}

##Julietella oli vaarin nama, tama funktio kaantaa ne oikein. VAARALLINEN MUTKA!!!
convertBackAlleleStatus <- function(x){
  if (x == 0){
    return('A A')
  }
  if (x == 1){
    return('A P')
  }
  if (x ==2){
    return('P P')
  }
  else{
    return(x)
  }
}

buildPresenceAbsenceByPosition <- function(positions, mapWithNames){
        #keep only first of the duplicated(at lower resolution)
        #HLA in the position file
        positions$P.group <- NULL
        positions_no_duplicates <- positions[!duplicated(positions$Trimmed), ]
        position_dictionary <- sapply(paste0(positions_no_duplicates$Locus, positions_no_duplicates$Trimmed), 
                                      function(x){gsub(':', '', x)}, USE.NAMES = FALSE)
        row.names(positions_no_duplicates) <- position_dictionary
        positions_no_duplicates <- positions_no_duplicates[, -c(1:4)]
        positions_no_duplicates <- positions_no_duplicates[which(row.names(positions_no_duplicates) %in% names(mapWithNames)), ] ##removing unnecessary alleles to simplify further computation
        positions_no_duplicates <- positions_no_duplicates[, apply(positions_no_duplicates, 2, function(x){!(length(unique(x)) == 1)})] ##keep only ones with variation
        long_names <- character()
        for (name in names(positions_no_duplicates)){
          long_names <- c(long_names, paste0(name, unique(unlist(positions_no_duplicates[name], use.names = FALSE))))
        }
        
        long_name_df <- data.frame(matrix(rep(0, length(long_names)*nrow(positions_no_duplicates)), nrow = nrow(positions_no_duplicates), ncol = length(long_names)))
        names(long_name_df) <- long_names
        row.names(long_name_df) <- row.names(positions_no_duplicates)
        
        for(column in names(positions_no_duplicates)){
          for(row in row.names(positions_no_duplicates)){
            long_name_df_colname <- paste0(column, positions_no_duplicates[row, column])
            long_name_df[row, long_name_df_colname] <- long_name_df[row, long_name_df_colname] + 1    
          }    
        }
        mapWithNames <- mapWithNames[, -c(1:6)]
        mapWithNames <- data.frame(apply(mapWithNames, c(1,2), convertAlleleStatus), stringsAsFactors = FALSE) 
        resultdf <- data.frame(matrix(vector(), 0, ncol(long_name_df), dimnames = list(c(), names(long_name_df))), stringsAsFactors = FALSE)
        
        for(rowindex in 1:nrow(mapWithNames)){
          nonzeronames <- names(mapWithNames)[which(mapWithNames[rowindex, ] != 0)]
          if(length(nonzeronames) == 1){
            result <- long_name_df[nonzeronames[[1]], ]*2
          }
          else{
            result <- long_name_df[nonzeronames[[1]], ] + long_name_df[nonzeronames[[2]], ]
          }
          resultdf <- rbind(resultdf, result)
        }
        row.names(resultdf) <- NULL
        return(resultdf)
                
}
##Special case below
# test <- buildMapWithNames(pedfile = "HLA-alleles-in-Ichip-cases.ped", mapfile = "HLA-alleles-in-Ichip-cases.map" )
# test <- test[, c(1:6,grep("^A", names(test)))] ##A only
# 
# 
# positions <- read.table("ExonPtnAlign_A.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = c(' ', '*', 'NA', '0 0'))
# positions$P.group <- NULL #!!!!!
# positions_no_duplicates <- positions[!duplicated(positions$Trimmed), ]
# posdict <- sapply(paste0(positions_no_duplicates$Locus, positions_no_duplicates$Trimmed), 
#                   function(x){gsub(':', '', x)}, USE.NAMES = FALSE)
# 
# row.names(positions_no_duplicates) <- posdict
# positions_no_duplicates <- positions_no_duplicates[, -c(1:4)]
# positions_no_duplicates <- positions_no_duplicates[which(row.names(positions_no_duplicates) %in% names(test)), ] ##removing unnecessary alleles to simplify further computation
# positions_no_duplicates <- positions_no_duplicates[, apply(positions_no_duplicates, 2, function(x){!(length(unique(x)) == 1)})] ##keep only ones with variation
# 
# long_names <- character()
# for (name in names(positions_no_duplicates)){
#   long_names <- c(long_names, paste0(name, unique(unlist(positions_no_duplicates[name], use.names = FALSE))))
# }
# 
# long_name_df <- data.frame(matrix(rep(0, length(long_names)*nrow(positions_no_duplicates)), nrow = nrow(positions_no_duplicates), ncol = length(long_names)))
# names(long_name_df) <- long_names
# row.names(long_name_df) <- row.names(positions_no_duplicates)
# 
# for(column in names(positions_no_duplicates)){
#   for(row in row.names(positions_no_duplicates)){
#     long_name_df_colname <- paste0(column, positions_no_duplicates[row, column])
#     long_name_df[row, long_name_df_colname] <- long_name_df[row, long_name_df_colname] + 1    
#   }    
# }
# 
# test <- test[, -c(1:6)]
# test <- data.frame(apply(test, c(1,2), convertAlleleStatus), stringsAsFactors = FALSE)
# 
# resultdf <- data.frame(matrix(vector(), 0, ncol(long_name_df), dimnames = list(c(), names(long_name_df))), stringsAsFactors = FALSE)
# 
# for(rowindex in 1:nrow(test)){
#   nonzeronames <- names(test)[which(test[rowindex, ] != 0)]
#   if(length(nonzeronames) == 1){
#     result <- long_name_df[nonzeronames[[1]], ]*2
#   }
#   else{
#     result <- long_name_df[nonzeronames[[1]], ] + long_name_df[nonzeronames[[2]], ]
#   }
#   resultdf <- rbind(resultdf, result)
# }
# 


###USAGE:
  # test <- buildMapWithNames(pedfile = "HLA-alleles-in-Ichip-cases.ped", mapfile = "HLA-alleles-in-Ichip-cases.map" )
  # test <- test[, c(1:6,grep("^A", names(test)))] ##A only
  # positions <- read.table("ExonPtnAlign_A.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = c(' ', '*', 'NA', '0 0'))
  #test <- data.frame(test,buildPresenceAbsenceByPosition(positions, test), stringsAsFactors = FALSE)