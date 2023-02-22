import sys
import yaml
import subprocess

def parse_config(config, encoding='utf-8'):
	with open(config) as file:
		yamldict = yaml.load(file, Loader=yaml.FullLoader)
	return yamldict
	
	
def run_associations(yamldict):
	associations = yamldict['assoc_list']
	digits       = yamldict['digits']
	for assoc in associations:
		print(assoc)
		inputfile  = yamldict['workfolder']+assoc+'/'+assoc+'.input'
		covarfile  = yamldict['workfolder']+assoc+'/'+assoc+'.covar'

		for model in yamldict['tests']:
			test        = yamldict['tests'][model]
			outputfile  = yamldict['workfolder']+assoc+'/'+assoc+'.'+model+'.'+test+'.assoc'
			summaryfile = yamldict['workfolder']+assoc+'/'+assoc+'.'+model+'.'+test+'.summary'
			# run summary statistics
			cmd        = ['python', 'PyHLA.py', '--input', inputfile,
			              '--summary', '--out', summaryfile,
						  '--digit', str(digits)]
			try:
				subprocess.check_call(cmd)
			except subprocess.CalledProcessError:
				print('Summary statistics failed on '+assoc+', '+model+', '+test)
				sys.exit(1)

			# run association analysis
			cmd        = ['python', 'PyHLA.py', '--input', inputfile,
			              '--assoc', '--model', model, '--test', test,
						  '--digit', str(digits), '--adjust', 'FDR',
						  '--perm', '100', '--seed', '42',
						  '--out', outputfile, '--covar', covarfile]
			try:
				subprocess.check_call(cmd)
			except subprocess.CalledProcessError:
				print('Association failed on '+assoc+', '+model+', '+test)
				sys.exit(1)

def main():
	yamldict = parse_config(sys.argv[1])
	run_associations(yamldict)


if __name__ == '__main__':
	main()
