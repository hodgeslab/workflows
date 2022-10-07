to make genome browser tracks using pyGenomeTracks, log in to the cluster via ssh, then run the following:

module load python/3.7.3
module load bedtools
pip3 --user install pyGenomeTracks

The system will install the pyGenomeTracks python package and its dependencies. Assuming installation completes successfully, you can run the program by running:

pyGenomeTracks

If you log in and log out, you will need to reload the python3 and bedtools environment modules each time before you can run it again:

module load python/3.7.3
module load bedtools

You will then need to edit the "tracks.ini" file (example attached) and link to it in your path. This sets up the details of the individual tracks. Edit the tracks.ini file to match what you want to create for your figure, such as the color, min and max limits, add more, remove things, etc. Then place it in your folder and run:

pyGenomeTracks --tracks tracks.ini \
                 --region chr1:1000000-4000000 \
                 --plotWidth 2 \
                 --out tracks.pdf \
                 --fontSize 7

Replace the "region" argument with the actual genome region you want to plot. Some parameters can also be adjusted, for example the plotWidth, to make the region wider or more narrow. Don't worry about the text labels in the output, these should be fixed in Illustrator.

You can find overall documentation here: https://pygenometracks.readthedocs.io and other examples for the tracks.ini file here: https://pygenometracks.readthedocs.io/en/latest/content/examples.html#basic-examples
