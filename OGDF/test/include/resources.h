/** \file
 * \brief Helper functions to be used in tests for accessing resource files.
 *
 * \author Tilo Wiedera
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_TEST_RESOURCES
#define OGDF_TEST_RESOURCES

#include <bandit/bandit.h>
#include <tinydir.h>
#include <vector>

#include <ogdf/fileformats/GraphIO.h>

using namespace bandit;
using namespace ogdf;
using std::string;
using std::function;
using std::vector;

const string RESOURCE_DIR = "test/resources";
typedef bool (*GraphReader)(Graph&, const char*);

/**
 * Tests whether the resource directory is present (i.e. the working directory is correct).
 *
 * \return true iff the resource directory was found
 */
inline bool resourceCheck() {
	tinydir_dir dir;
	bool result = tinydir_open(&dir, RESOURCE_DIR.c_str())  != -1;
	tinydir_close(&dir);

	return result;
}
/**
 * Iterates over each file contained in the specified directory.
 *
 * \param directory The path of the directory.
 * \param callback A function that will be called for each file in the directory.
 * \param recurse Whether to include sub directories.
 */
inline void for_each_file(const string &directory, function<void (const string&)> callback, bool recurse = false) {
	string resourceDirectory = RESOURCE_DIR + "/" + directory;
	tinydir_dir dir;

	if(tinydir_open(&dir, resourceDirectory.c_str())  == -1) {
	    	it("", [&](){ Assert::Failure("Could not open directory: " + resourceDirectory); });
	} else while (dir.has_next) {
	    tinydir_file file;

	    if(tinydir_readfile(&dir, &file) == -1) {
	    	it("", [&](){ Assert::Failure("Could not read directory: " + resourceDirectory); });
	    	break;
	    }

	    if (!file.is_dir) {
	    	callback(resourceDirectory + "/" + file.name);
	    } else if(recurse) {
	    	for_each_file(directory + "/" + file.name, callback, true);
	    }

	    tinydir_next(&dir);
	}

	tinydir_close(&dir);
}

/**
 * Reads the specified files and creates a test for each graph.
 *
 * \param title The base title for the test cases.
 * \param callback filenames The names of the files to be read.
 * \param recurse testFunc The actual test to be performed.
 * \param reader The function used to parse the files, defaults to GraphIO::readGML.
 */
inline void for_each_graph_it(const string &title, const vector<string> &filenames, function<void (Graph &graph, const string&)> testFunc, GraphReader reader = GraphIO::readGML) {
	for(const string filename : filenames) {
		it(string(title + " [" + filename.c_str() + "] ").c_str(), [&](){
			Graph graph;
			AssertThat(reader(graph, (RESOURCE_DIR + "/" + filename).c_str()), IsTrue());
			testFunc(graph, filename);
		});
	}
}

#endif
