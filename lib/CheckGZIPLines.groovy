#!/usr/bin/env groovy
import java.util.zip.GZIPInputStream
import java.io.IOException

static def checkFastqFile(filePath) {
    File file = new File(filePath)

    if (!file.exists()) {
        println "File not found: ${filePath}"
        return false
    }

    try {
        GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(file))
        BufferedReader reader = new BufferedReader(new InputStreamReader(gzipInputStream))

        int lineCount = 0

        while ((reader.readLine()) != null) {
            lineCount++
            if (lineCount >= 4) {
                reader.close()
                return true
            }
        }
        reader.close()
    } catch (IOException e) {
        println "Error reading file: ${e.message}"
        return false
    }

    return false
}

static def processMeta( filePaths ) {
    if (!(filePaths instanceof List)) {
        filePaths = filePaths.toList()
    }

    boolean pass = filePaths.every { filePath ->
        checkFastqFile(filePath as String)
    }

    return [filePaths, pass]
}