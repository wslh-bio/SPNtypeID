#!/usr/bin/env groovy

import java.nio.file.Files
import java.nio.file.Paths

// // Function to check if the first four lines of a gzipped file have content
// boolean checkFileContent(String filePath) {
//     def process = ["zcat", filePath].execute() // Use zcat to read gzipped content
//     def content = process.in.withReader { it.readLines() } // Read all lines from the stream
    
//     // Ensure the first four lines have content
//     return content.take(4).any { it?.trim() }
// }

// // Main function to process samples
// Map processSamples(Map<String, List<String>> samples, String logFilePath) {
//     def logFile = new File(logFilePath)
//     logFile.text = "" // Clear any existing content

//     Map validSamples = [:]

//     samples.each { sampleName, fastqFiles ->
//         def isValid = fastqFiles.every { file ->
//             if (!checkFileContent(file)) {
//                 logFile << "File is empty: ${file}\n"
//                 return false // Skip this sample if any file is empty
//             }
//             return true
//         }
//     }  

class FileStructures (MetaID, Paths) {

    public String MetaID
    public String Paths
    public String Empty

    def setMetaID(String MetaID) {
        MetaID = name
    }

    def setPaths(String Paths) {
        Paths = files
    }

    def setContents(String Empty) {
        Empty = content
    }

    def getFile() {
        println("$paths")
    }
    
    def getName() {
        println("$MetaID")
    }

    def checkContents(Paths)

    public static void main(args) {
        FileStructures =  setMetaID(MetaID, Paths)
    }

}