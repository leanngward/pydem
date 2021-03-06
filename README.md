# pydem
# Creating PyDem Virtual Environment

## Initial Set Up Requirements:

1. Python >= 3.7 <br>
2. pip installer <br>
3. virtualenv <br>


To install the virtualenv, use the following line of code in the command window: <br>
`$ pip install --user virtualenv`
	                                                                  

# Running PyDem:

1. Create a new environment using `$ venv myenvname`
2. Activate environment using `.\myenvname\Scripts\activate`
3. Save requirements.txt files to desired folder. 
4. Install needed environment using ` $ pip install -r requirements.txt` 
5. Run desired project file inside virutal environment. (Note: Be sure to have pydem.py and the project file saved in the same location. Set your path to this location in the command line or move pydem.py to the python path if applicable.)


#### Future Notes:
1. New packages that are needed can be installed with your environment activated using pip install or pip --upgrade. 
2. If environment is updated in the future, the requirements.txt file can be exported again using `pip freeze` while the environment is activated. 
This requirements.txt file can be shared using the list of steps above.
