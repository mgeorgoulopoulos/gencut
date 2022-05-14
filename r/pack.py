import re

files = ['includes.cpp'
	, '../src/utils/AtomBox.h' 
	, '../src/utils/Exception.h' 
	, '../src/utils/Vec3D.h' 
	, '../src/app/StringConstants.cpp' 
	, '../src/app/GeneCollections.h' 
	, '../src/app/GeneSignal.h' 
	, '../src/app/GeneSignal.cpp' 
	, '../src/app/GeneSignalMatrix.h' 
	, '../src/app/GeneSignalMatrix.cpp' 
	, '../src/app/GenomeModel.h' 
	, '../src/app/GenomeModel.cpp' 
	, '../src/app/GenomeCutter.h' 
	, '../src/app/GenomeCutter.cpp' 
	, 'rapi.cpp'
	]

all = '';

for file in files:
	# Tell filename
	all = all + '// ---------------------------------\n'
	all = all + '// Start file: ' + file + '\n'
	all = all + '// ---------------------------------\n'
	all = all + '\n'
	
	# Load and append file
	f = open(file,mode='r')
	data = f.read()
	f.close()
	
	# Remove headers
	if file != 'includes.cpp':
		data = re.sub('#include.*\n', '', data)
	
	all += data
	all += '\n'




# Write big file
f = open('gencut-packed.cpp',mode='w')
f.write(all)
f.close()
