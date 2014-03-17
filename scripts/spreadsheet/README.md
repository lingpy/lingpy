Spreadsheet reader examples
===========================

This directory contains data and an example script to load CSV data (output from a spreadsheet) into a LingPy wordlist object and to output the data as QLC format.

The data come from Heath et al 2014: [http://dogonlanguages.org/lexicons.cfm] (http://dogonlanguages.org/lexicons.cfm)

- Heath2014.prf - orthography profile that specifies the Unicode grapheme clusters in the Dogon wordlist
- dogon.bl - a blacklist file that contains triple-comma delimited <,,,> find and replace regular expressions
- dogon_wordlists.tsv - a subset of the Heath et al 2014 Dogon Comparative Wordlist
- heath2014-wordlist.py - script illustrating how to load the Dogon comparative wordlist and write QLC data format


@misc{Heath_etal2014,
	Author = {Jeffrey Heath and Steven Moran and Laura McPherson and Kirill Prokhorov and Abbie Hantgan and Brian Cansler},
	Howpublished = {Online: \url{http://dogonlanguages.org}},
	Title = {Dogon Comparative Wordlist},
	Year = {2014}}

