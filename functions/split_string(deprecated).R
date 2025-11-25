#takes a string like JS0002-20110427-001-1.wav
#and splits it into recording and selec. 
#i.e. JS0002-20110427-001.wav and 1

split_string <- function(string) {
  # Split the string at the last '-'
  parts <- strsplit(string, "-")[[1]]
  # Extract the last part
  last_part <- tail(parts, 1)
  # Split the last part into the string and the number
  split_last_part <- strsplit(last_part, "\\.")[[1]]
  extracted_string <- paste(parts[-length(parts)], collapse = "-")
  #add .wav to match sound.files syntax
  extracted_string <- paste(extracted_string, ".wav", sep="")
  extracted_number <- as.numeric(split_last_part[1])
  return(list( sound.files= extracted_string, selec = extracted_number))
}

# Test the function
string <- "JS0002-20110427-001-1.wav"
result <- split_string(string)
print(result$string)
print(result$number)

#add proper tests