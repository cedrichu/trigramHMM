import operator
import sys
import optparse

def replace_RARE(output_file):
	train_file = open("gene.train","r")
	output_file = open(output_file, "w")
	
	word_counter = count_words()

	for line in train_file:
		if line == "\n":
			output_file.write("\n")
			continue
		line_buffer = line.split()
		if word_counter[line_buffer[0]] < 5:
			output_file.write("_RARE_ "+line_buffer[1]+"\n")
		else:
			output_file.write(line)

def replace_RARE_improvement(output_file):
	train_file = open("gene.train","r")
	output_file = open(output_file, "w")
	
	word_counter = count_words()

	for line in train_file:
		if line == "\n":
			output_file.write("\n")
			continue
		line_buffer = line.split()
		if word_counter[line_buffer[0]] < 5:

			if check_numeric(line_buffer[0]):
				output_file.write("_NUMERIC_ "+line_buffer[1]+"\n")
			elif check_all_capitals(line_buffer[0]):
				output_file.write("_ALL_CAPITALS_ "+line_buffer[1]+"\n")
			elif check_last_capital(line_buffer[0]):
				output_file.write("_LAST_CAPITAL_ "+line_buffer[1]+"\n")
			else:
				output_file.write("_RARE_ "+line_buffer[1]+"\n")
		else:
			output_file.write(line)


def count_words():
	count_file = open("gene.counts","r")
	word_counter = {}
	for line in count_file:
		line_buffer = line.split()
		if line_buffer[1] == "WORDTAG":
			if line_buffer[3] in word_counter:
				word_counter[line_buffer[3]] += int(line_buffer[0])
			else:
				word_counter[line_buffer[3]] = int(line_buffer[0])
	return word_counter

def check_numeric(word):
	return any(char.isdigit() for char in word)

def check_all_capitals(word):
	return all(char.isupper() for char in word)

def check_last_capital(word):
	return any(char.islower() for char in word[:-1]) and word[-1].isupper()


def get_options(args=None):
    optParser = optparse.OptionParser()
    optParser.add_option("-t", "--train-file", dest="train_file",
                            help="define the output train file (mandatory)")
    optParser.add_option("-i", action="store_true",
                         default=False, help="if add improvement in problem3")
    (options, args) = optParser.parse_args(args=args)
    return options

if __name__ == "__main__":
	options = get_options()
	if options.i == False:
		replace_RARE(options.train_file)
	else:
		replace_RARE_improvement(options.train_file)





