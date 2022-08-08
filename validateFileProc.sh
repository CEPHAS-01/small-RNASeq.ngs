#!/bin/bash
function uniValidateFileProc()
{
        #Use the function 'uniValidateFile and produce loggable output for the log file
        outCode=`uniValidateFile $1 $2`
#        echo $outCode
        if [[ $outCode = 10 ]];
        then
                echo "yes"

        elif [[ $outCode = 5 ]];
        then
                echo "File $1 does not exist!"

        elif [[ $outCode = 6 ]];
        then
                echo "The extension of the file $1 is not appropriate"

        elif [[ $outCode = 7 ]];
        then
                echo "File $1 is not readable"

        elif [[ $outCode = 8 ]];
        then
                echo "File $1 is empty"

        elif [[ $outCode = 9 ]];
        then
                echo "This is not a directory"

	elif [[ $outCode = "incomplete" ]];
        then
                echo "Incorrect number of parameters entered for function 'uniValidateFileProc' !"

        else
                echo "nil"
        fi
}
