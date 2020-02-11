@ECHO OFF
SETLOCAL
IF EXIST D:\wamp\www\pm1nps\courses\MAS201 (
 SET WEBDIR=D:\wamp\www\pm1nps\courses\MAS201
) ELSE (
IF EXIST C:\wamp\www\courses\MAS201 (
 SET WEBDIR=C:\wamp\www\courses\MAS201
) ELSE (
 SET WEBDIR=C:\wamp\www\pm1nps\courses\MAS201
))

echo "WEBDIR="
echo %WEBDIR%

pdflatex all_lectures
pdflatex -job-name=all_handouts \def\HO{1} \input{all_lectures}

copy all_lectures.pdf %WEBDIR%\lectures\all_lectures.pdf
copy all_handouts.pdf %WEBDIR%\lectures\all_handouts.pdf

