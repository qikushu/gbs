## Usage
# Download these functions from the GitHub repository
# source("https://raw.githubusercontent.com/qikushu/gbs/master/drawGraphicalGenotype.R")

# Read data (do not run)
# plantGenoList = readGGdata("test.csv")

# Plot and save in PNG file (do not run)
# png("output.png", width = 2400, height = 2000)
# plotHorizontal(plantGenoList = plantGenoList, chrMaxLength = chrMaxLength, color_mapping = color_mapping)
# dev.off()

## Other parameters
# Total rice chromosome length of Nipponbare assembly in Mb
chrMaxLength <- c("1" = 43.270923, "2" = 35.937250, "3" = 36.413819, "4" = 35.502694, "5" = 29.958434, "6" = 31.248787, "7" = 29.697621, "8" = 28.443022, "9" = 23.012720, "10" = 23.207287, "11" = 29.02110612, "12" = 27.531856)

# Conversion rule of ABH genotypes and color code in output graphics.
color_mapping <- c("A" = "orange", "B" = "#00bfff", "H" = "red", "-" = "gray")

## Functions
readGGdata = function(infile) {

    # Read data
    rawdf = read.csv(infile, header = TRUE)

    # Get plant names (Plant column)
    plantNames = rawdf$Plant[-c(1, 2)]

    # Evaluate the first row to get non-empty column indices
    notEmptyIndex = which(rawdf[1, ] != "")
    df1 = rawdf[, notEmptyIndex]
    genotypeData = as.data.frame(t(df1))
    colnames(genotypeData) = append(c("Chr", "Pos"), plantNames)
    genotypeData$Pos = as.numeric(genotypeData$Pos)  # Convert to numeric type
    genotypeData$Marker = rownames(genotypeData)

    # Split data for each plant
    forDeletion = c("Marker", "Chr", "Pos")
    plantGenoList <- lapply(plantNames, function(col_name) {
        cols_to_include <- c(forDeletion, col_name)
        df_subset <- genotypeData[cols_to_include]
        colnames(df_subset) <- cols_to_include
        return(df_subset)
    })
    # Hold data for each plant in a named list
    names(plantGenoList) <- plantNames

    return(plantGenoList)
}

plotHorizontal = function(plantGenoList, chrMaxLength, color_mapping) {

    ### Code for arranging chromosomes horizontally
    # Spacing between chromosomes
    chrInterval = 5
    # Spacing between plants
    plantInterval = 2
    # Width of horizontal bars for markers
    width = 12

    # Left margin spacing
    leftSpace = 10
    # Right margin spacing
    rightSpace = 10
    # Top margin spacing
    topSpace = 0

    # Hold plant names in a named list
    plantNames = names(plantGenoList)
    # Number of plants
    plantNum = length(plantNames)
    # Chromosome names
    chrNames = unique(plantGenoList[[1]]$Chr)
    # Number of chromosomes
    numChr = length(chrNames)
    # All chromosome names
    allChrNames = names(chrMaxLength)
    # Number of all chromosomes
    allChrNum = length(allChrNames)

    # Maximum width for the x-axis
    maxXwidth = leftSpace + sum(chrMaxLength) + numChr * chrInterval + rightSpace + 10

    # Maximum height for the y-axis
    maxYwidth = (width + plantInterval) * plantNum

    # Set up the plot area
    plot(1, type = "n", xlim = c(-20, maxXwidth), ylim = c(0, maxYwidth), xlab = "X-axis", ylab = "Y-axis")

    # Process each plant's data
    j = 1
    # Offset in the y-direction for chromosome drawing
    yoff = 0;

    # Loop over the named list `chromosome_data`
    for (plantName in plantNames) {
        cat(plantName)
        # Get data for each plant
        df <- plantGenoList[[plantName]]
        colnames(df)[colnames(df) == plantName] <- "Genotype"
        # Offset in the x-direction for chromosome drawing
        xoff = 0;
        # Calculate y_start coordinates
        y_start = topSpace + yoff
        # Calculate y_end coordinates
        y_end = topSpace + yoff + width
        # Calculate y_mid
        y_mid = (y_start + y_end) / 2
        # Write the lineage name
        text(x = -20, y = y_start, labels = plantName, cex = 1.5)

        for (chr in allChrNames) {

            # End of chromosome
            chrEnd = chrMaxLength[chr]

            # Check if user data exists
            if (any(df$Chr == chr)) {
                # Extract data for each chromosome
                df2 = df[df$Chr == chr, ]

                # Handling when adjacent genotypes are the same
                df3 = flankGeno(df2, chr)

                # Add a new column "color" to store color information
                df3$color <- color_mapping[df3$Genotype]

                # Loop through the data frame df3
                for (i in 1:nrow(df3)) {
                    row <- df3[i, ]  # Extract one row
                    # Array of xleft coordinates
                    xleft <- df3$MinPos + xoff
                    # Array of xright coordinates
                    xright <- df3$MaxPos + xoff
                    # Genotype color
                    color = df3$color
                    rect(xleft = xleft, ytop = y_start, xright = xright, ybottom = y_end, col = color, border = 0)
                }
            }
            rect(xleft = 0 + xoff, ytop = y_start, xright = chrEnd + xoff, ybottom = y_end, col = NA, border = "grey")
            xoff = xoff + chrEnd + chrInterval
        }
        yoff = yoff + width + plantInterval
        j = j + 1;
    }
}

# Function to return a data frame with adjacent genotypes grouped together
flankGeno = function(df2, chr) {
    df2$SequenceNumber <- createSequenceDataFrame(df2$Genotype)$SequenceNumber
    # Maximum value of SequenceNumber
    maxLinkGroup = max(df2$SequenceNumber)

    result_df <- data.frame()  # Initialize an empty data frame

    for (i in 1:maxLinkGroup) {
        df3 = df2[df2$SequenceNumber == i, ]
        # Genotype
        gt = df3$Genotype[1]
        maxpos = max(df3$Pos)
        minpos = min(df3$Pos)
        result_df <- rbind(result_df, data.frame(Chr = chr, MinPos = minpos, MaxPos = maxpos, Genotype = gt))
    }

    return(result_df)
}

# Assign sequential numbers to adjacent genotypes when they are the same
createSequenceDataFrame <- function(alphabet_vector) {
    # Initialize a vector to store sequential numbers
    sequence_number <- rep(NA, length(alphabet_vector))

    # Set the first element to sequence number 1
    sequence_number[1] <- 1

    # Assign sequential numbers
    for (i in 2:length(alphabet_vector)) {
        if (alphabet_vector[i] == alphabet_vector[i - 1]) {
            sequence_number[i] <- sequence_number[i - 1]
        } else {
            sequence_number[i] <- sequence_number[i - 1] + 1
        }
    }

    # Return the result as a data frame
    result_df <- data.frame(Alphabet = alphabet_vector, SequenceNumber = sequence_number)
    return(result_df)
}
