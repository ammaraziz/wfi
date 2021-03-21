#!/bin/bash -x

if [ -r ../LABEL ];then
	LABEL=../LABEL
else
	LABEL=$(which LABEL)
fi


if [ -r ../IRMA ];then
	IRMA=../IRMA
else
	IRMA=$(which IRMA)
fi

$LABEL test1.fa label-test H9v2011
$IRMA FLU test2.fastq irma-test
