#ifndef XMLTOOLS_H
#define XMLTOOLS_H
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/dom/DOM.hpp>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/AbstractDOMParser.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSParser.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMError.hpp>
#include <xercesc/dom/DOMLocator.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMAttr.hpp>
#include <Eigen/Core>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <iostream>

USING_PART_OF_NAMESPACE_EIGEN

XERCES_CPP_NAMESPACE_USE

// ---------------------------------------------------------------------------
//  Simple error handler deriviative to install on parser
// ---------------------------------------------------------------------------
class DOMCustomErrorHandler : public DOMErrorHandler
{
public:
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    DOMCustomErrorHandler();
    ~DOMCustomErrorHandler();


    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    bool getSawErrors() const;


    // -----------------------------------------------------------------------
    //  Implementation of the DOM ErrorHandler interface
    // -----------------------------------------------------------------------
    bool handleError(const DOMError& domError);
    void resetErrors();


private :
    // -----------------------------------------------------------------------
    //  Unimplemented constructors and operators
    // -----------------------------------------------------------------------
    DOMCustomErrorHandler(const DOMCustomErrorHandler&);
    void operator=(const DOMCustomErrorHandler&);


    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fSawErrors
    //      This is set if we get any errors, and is queryable via a getter
    //      method. Its used by the main code to suppress output if there are
    //      errors.
    // -----------------------------------------------------------------------
    bool    fSawErrors;
};


// ---------------------------------------------------------------------------
//  This is a simple class that lets us do easy (though not terribly efficient)
//  trancoding of XMLCh data to local code page for display.
// ---------------------------------------------------------------------------
class StrX
{
public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    StrX(const XMLCh* const toTranscode)
    {
        // Call the private transcoding method
        fLocalForm = XMLString::transcode(toTranscode);
    }

    ~StrX()
    {
        XMLString::release(&fLocalForm);
    }


    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const char* localForm() const
    {
        return fLocalForm;
    }

private :
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fLocalForm
    //      This is the local code page form of the string.
    // -----------------------------------------------------------------------
    char*   fLocalForm;
};

inline XERCES_STD_QUALIFIER ostream& operator<<(XERCES_STD_QUALIFIER ostream& target, const StrX& toDump)
{
    target << toDump.localForm();
    return target;
}

inline bool DOMCustomErrorHandler::getSawErrors() const
{
    return fSawErrors;
}

class XMLTools
{
private:
	XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *m_doc;
	DOMLSParser       *m_parser;
public:
	XMLTools();
	~XMLTools();
	bool loadFile(std::string filename);

	bool createFile(std::string masterTag);
	bool saveToFile(std::string filename);
	void closeHandle() { m_doc->release(); }
	bool handleData(bool write, std::string &data, std::string tagName, int itemId = 0);
	bool handleData(int &data, const std::string &tagName, const std::string &attributeName, const int &itemId = 0);
	bool handleData(float &data, const std::string &tagName, const std::string &attributeName, const int &itemId = 0);
	bool handleData(bool write, std::string &data, const std::string &tagName, const std::string &attributeName, const int &itemId = 0);
	int numberOfElements(std::string tagName);
	bool handleData(Eigen::Matrix4f &matrix, std::string tagName, int itemId = 0);

	/**
	@brief Helper function that parses a string an outputs an Eigen Library matrix
	**/
	void string2matrix(Eigen::Matrix4f &matrix, std::string &str);
};

#endif