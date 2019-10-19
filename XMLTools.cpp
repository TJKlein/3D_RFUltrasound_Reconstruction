#include "XMLTools.h"
#include <boost/tokenizer.hpp> 

/* Note:
Problems with unresolved external: Make sure under project properties that Treat wchar_t as
Built-in Type is set to Yes (C/C++, the language tab).
*/

#include <QtWidgets/QMessageBox>
#include <sstream>
XMLTools::XMLTools() :
m_doc(0),
m_parser(0)
{}

XMLTools::~XMLTools()
{

	//
	//  Delete the parser itself.  Must be done prior to calling Terminate, below.
	//

	if ( m_parser)
		m_parser->release();

	// And call the termination method
	XMLPlatformUtils::Terminate();
}



bool XMLTools::saveToFile(std::string xmlFile)
{

	bool                       errorOccurred = false;

	try
	{
		AbstractDOMParser::ValSchemes valScheme = AbstractDOMParser::Val_Never;

		// initialize if not previously done
		if ( m_parser == 0 )
			XMLPlatformUtils::Initialize();
	}

	catch (const XMLException& toCatch)
	{
		XERCES_STD_QUALIFIER cerr << "Error during initialization! :\n"
			<< StrX(toCatch.getMessage()) << XERCES_STD_QUALIFIER endl;
		return false;
	}
	// Instantiate the DOM parser.
	static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
	DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(gLS);
	// check if not previously initialized
	if ( m_parser == 0)
		m_parser = ((DOMImplementationLS*)impl)->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
	DOMConfiguration  *config = m_parser->getDomConfig();
	DOMLSSerializer   *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer();
    DOMLSOutput       *theOutputDesc = ((DOMImplementationLS*)impl)->createLSOutput();


	config->setParameter(XMLUni::fgDOMNamespaces, false);
	config->setParameter(XMLUni::fgXercesSchema, false);
	config->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
	config->setParameter(XMLUni::fgDOMValidate, false);

	// enable datatype normalization - default is off
	config->setParameter(XMLUni::fgDOMDatatypeNormalization, false);

	// And create our error handler and install it
	DOMCustomErrorHandler errorHandler;
	config->setParameter(XMLUni::fgDOMErrorHandler, &errorHandler);


	unsigned long duration;


	//reset error count first
	errorHandler.resetErrors();



	try
	{
		// plug in user's own error handler
            DOMConfiguration* serializerConfig=theSerializer->getDomConfig();

			XMLFormatTarget *myFormTarget = new LocalFileFormatTarget(xmlFile.c_str());

			theOutputDesc->setByteStream(myFormTarget);

			  theSerializer->write(m_doc, theOutputDesc);

            theOutputDesc->release();
            theSerializer->release();

	}
	catch (const XMLException& toCatch)
	{
		XERCES_STD_QUALIFIER cerr << "\nError writing: '" << xmlFile.c_str() << "'\n"
			<< "Exception message is:  \n"
			<< StrX(toCatch.getMessage()) << "\n" << XERCES_STD_QUALIFIER endl;
		errorOccurred = true;
	}
	catch (const DOMException& toCatch)
	{
		const unsigned int maxChars = 2047;
		XMLCh errText[maxChars + 1];

		XERCES_STD_QUALIFIER cerr << "\nDOM Error during writing: '" << xmlFile.c_str() << "'\n"
			<< "DOMException code is:  " << toCatch.code << XERCES_STD_QUALIFIER endl;

		if (DOMImplementation::loadDOMExceptionMsg(toCatch.code, errText, maxChars))
			XERCES_STD_QUALIFIER cerr << "Message is: " << StrX(errText) << XERCES_STD_QUALIFIER endl;

		errorOccurred = true;
	}
	catch (...)
	{
		XERCES_STD_QUALIFIER cerr << "\nUnexpected exception during writing: '" << xmlFile.c_str() << "'\n";
		errorOccurred = true;
	}


	return true;
}




bool XMLTools::loadFile(std::string xmlFile)
{
	bool                       errorOccurred = false;

	try
	{
		AbstractDOMParser::ValSchemes valScheme = AbstractDOMParser::Val_Never;

		// initialize if not previously done
		if ( m_parser == 0 )
			XMLPlatformUtils::Initialize();
	}

	catch (const XMLException& toCatch)
	{
		XERCES_STD_QUALIFIER cerr << "Error during initialization! :\n"
			<< StrX(toCatch.getMessage()) << XERCES_STD_QUALIFIER endl;
		return false;
	}
	// Instantiate the DOM parser.
	static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
	DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(gLS);
	// check if not previously initialized
	if ( m_parser == 0)
		m_parser = ((DOMImplementationLS*)impl)->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
	DOMConfiguration  *config = m_parser->getDomConfig();


	config->setParameter(XMLUni::fgDOMNamespaces, false);
	config->setParameter(XMLUni::fgXercesSchema, false);
	config->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
	config->setParameter(XMLUni::fgDOMValidate, false);

	// enable datatype normalization - default is off
	config->setParameter(XMLUni::fgDOMDatatypeNormalization, false);

	// And create our error handler and install it
	DOMCustomErrorHandler errorHandler;
	config->setParameter(XMLUni::fgDOMErrorHandler, &errorHandler);


	unsigned long duration;


	//reset error count first
	errorHandler.resetErrors();



	try
	{
		// reset document pool
		m_parser->resetDocumentPool();

		const unsigned long startMillis = XMLPlatformUtils::getCurrentMillis();
		m_doc = m_parser->parseURI(xmlFile.c_str());
		const unsigned long endMillis = XMLPlatformUtils::getCurrentMillis();
		duration = endMillis - startMillis;
	}

	catch (const XMLException& toCatch)
	{
		XERCES_STD_QUALIFIER cerr << "\nError during parsing: '" << xmlFile.c_str() << "'\n"
			<< "Exception message is:  \n"
			<< StrX(toCatch.getMessage()) << "\n" << XERCES_STD_QUALIFIER endl;
		errorOccurred = true;
	}
	catch (const DOMException& toCatch)
	{
		const unsigned int maxChars = 2047;
		XMLCh errText[maxChars + 1];

		XERCES_STD_QUALIFIER cerr << "\nDOM Error during parsing: '" << xmlFile.c_str() << "'\n"
			<< "DOMException code is:  " << toCatch.code << XERCES_STD_QUALIFIER endl;

		if (DOMImplementation::loadDOMExceptionMsg(toCatch.code, errText, maxChars))
			XERCES_STD_QUALIFIER cerr << "Message is: " << StrX(errText) << XERCES_STD_QUALIFIER endl;

		errorOccurred = true;
	}
	catch (...)
	{
		XERCES_STD_QUALIFIER cerr << "\nUnexpected exception during parsing: '" << xmlFile.c_str() << "'\n";
		errorOccurred = true;
	}


	return true;
}

bool XMLTools::handleData(bool write, std::string &data, std::string tagName, int itemId)
{
	DOMNodeList* _elemList = m_doc->getDocumentElement()->getElementsByTagName(XMLString::transcode(tagName.c_str()));


	if ( _elemList->getLength() > 0 )
	{
		try
		{
			if ( !write)
				data = XMLString::transcode(_elemList->item(itemId)->getTextContent()); 
			else // read
				_elemList->item(itemId)->setTextContent(XMLString::transcode(data.c_str()));
		}
		catch(...)
		{
			return false;
		}
		return true; 
	}
	else // create the new tag 	
	{
		try
		{
			DOMElement* _rootElem = m_doc->getDocumentElement();
			DOMElement*  _newElem = m_doc->createElement(XMLString::transcode(tagName.c_str()));
			_newElem->setTextContent(XMLString::transcode(data.c_str()));
			_rootElem->appendChild(_newElem);
		}
		catch(...)
		{
			return false;
		}
		return true;
	}
}

int XMLTools::numberOfElements(std::string tagName)
{
	DOMNodeList* _elemList = m_doc->getDocumentElement()->getElementsByTagName(XMLString::transcode(tagName.c_str()));

	return _elemList->getLength();
}


bool XMLTools::handleData(float &data, const std::string &tagName, const std::string &attributeName, const int &itemId)
{

	DOMNodeList* _elemList = m_doc->getDocumentElement()->getElementsByTagName(XMLString::transcode(tagName.c_str()));

	if ( _elemList->getLength() > 0 )
	{
		try
		{
			if ( _elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str())))
			{
			std::string s = XMLString::transcode(_elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str()))->getTextContent());
			data = atof(s.c_str());
			}
			else
				return false;
		}
		catch(...)
		{
			return false;
		}
		return true; 
	}
	else
		return false;
}

bool XMLTools::handleData(int &data, const std::string &tagName, const std::string &attributeName, const int &itemId)
{

	DOMNodeList* _elemList = m_doc->getDocumentElement()->getElementsByTagName(XMLString::transcode(tagName.c_str()));

	if ( _elemList->getLength() > 0 )
	{
		try
		{
			std::string s = XMLString::transcode(_elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str()))->getTextContent());
			data = atoi(s.c_str());
		}
		catch(...)
		{
			return false;
		}
		return true; 
	}
	else
		return false;
}

bool XMLTools::handleData(bool write, std::string &data, const std::string &tagName, const std::string &attributeName, const int &itemId)
{

	DOMNodeList* _elemList = m_doc->getDocumentElement()->getElementsByTagName(XMLString::transcode(tagName.c_str()));
	std::string _tmp = "";
	if ( _elemList->getLength() > 0 && _elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str())) != NULL )
		_tmp = XMLString::transcode(_elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str()))->getTextContent());
	if ( (!write && _elemList->getLength() > 0 && _elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str())) != NULL) || (_elemList->getLength() > 0 && _elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str())) != NULL && _tmp.compare(data) == 0))
	{
		try
		{
			data = XMLString::transcode(_elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str()))->getTextContent());
		}
		catch(...)
		{
			return false;
		}
		return true; 
	}
	else // create the new tag 	
	{
		try
		{
			
			int _val = _elemList->getLength();
			bool _found = false;
			// check if the tag already exists ( if not  create a new node)
			/*for(int i=0;i<_elemList->getLength();i++)
			{
				
				if ( _elemList->item(i)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str())) != NULL )
				{
					std::string _tmp = XMLString::transcode(_elemList->item(itemId)->getAttributes()->getNamedItem(XMLString::transcode(attributeName.c_str()))->getTextContent());
					if (_tmp.compare(attributeName) == 0)
					{
						_found = true;
					}
				}

			}*/

			_found = (_elemList->getLength() > 0 );

			if ( !_found)
			{
			DOMElement* _rootElem = m_doc->getDocumentElement();
			DOMElement*  _newElem = m_doc->createElement(XMLString::transcode(tagName.c_str()));
			_newElem->setAttribute(XMLString::transcode(attributeName.c_str()), XMLString::transcode(data.c_str()));
			_rootElem->appendChild(_newElem);
			}
			else
			{
				DOMNode *_node = (_elemList->item(itemId));
				
				DOMElement *_element = dynamic_cast<DOMElement*>(_node);
			

				_element->setAttribute(XMLString::transcode(attributeName.c_str()), XMLString::transcode(data.c_str()));
			}
		}
		catch(...)
		{
			return false;
		}
		return true;
	}
}


void XMLTools::string2matrix(Eigen::Matrix4f &matrix, std::string &str)
{
	Eigen::Matrix4f _tempMatrix;
	boost::char_separator<char> sep(" \n\t"); 
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
			tokenizer tok(str, sep);
			float data[16];
			int counter = 0;
			for(tokenizer::iterator it = tok.begin();it != tok.end();++it)
			{
				data[counter++]=atof((*it).c_str());
			}

			_tempMatrix << Eigen::Map<Eigen::Matrix4f>(data);
			matrix = _tempMatrix.transpose();
}

bool XMLTools::handleData(Eigen::Matrix4f &matrix, std::string tagName, int itemId)
{
	Eigen::Matrix4f _tempMatrix;
	DOMNodeList* _elemList = m_doc->getDocumentElement()->getElementsByTagName(XMLString::transcode(tagName.c_str()));


	if ( _elemList->getLength() > 0 )
	{

		try
		{
			typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
			std::string s = XMLString::transcode(_elemList->item(itemId)->getTextContent()); 
			boost::char_separator<char> sep(" \n\t"); 
			tokenizer tok(s, sep);
			float data[16];
			int counter = 0;
			for(tokenizer::iterator it = tok.begin();it != tok.end();++it)
			{
				data[counter++]=atof((*it).c_str());
			}

			_tempMatrix << Eigen::Map<Eigen::Matrix4f>(data);
			matrix = _tempMatrix.transpose();
		}
		catch(...)
		{
			return false;
		}
		return true;  
	}
	else
		return false;
}

DOMCustomErrorHandler::DOMCustomErrorHandler() :

fSawErrors(false)
{
}

DOMCustomErrorHandler::~DOMCustomErrorHandler()
{
}


// ---------------------------------------------------------------------------
//  DOMCountHandlers: Overrides of the DOM ErrorHandler interface
// ---------------------------------------------------------------------------
bool DOMCustomErrorHandler::handleError(const DOMError& domError)
{
	fSawErrors = true;
	if (domError.getSeverity() == DOMError::DOM_SEVERITY_WARNING)
		XERCES_STD_QUALIFIER cerr << "\nWarning at file ";
	else if (domError.getSeverity() == DOMError::DOM_SEVERITY_ERROR)
		XERCES_STD_QUALIFIER cerr << "\nError at file ";
	else
		XERCES_STD_QUALIFIER cerr << "\nFatal Error at file ";

	std::stringstream ss;

	ss << StrX(domError.getLocation()->getURI())
		<< ", line " << domError.getLocation()->getLineNumber()
		<< ", char " << domError.getLocation()->getColumnNumber()
		<< "\n  Message: " << StrX(domError.getMessage()) << XERCES_STD_QUALIFIER endl;

	QMessageBox msgBox;
	msgBox.setText(ss.str().c_str());
	msgBox.exec();
	return true;
}

void DOMCustomErrorHandler::resetErrors()
{
	fSawErrors = false;
}


bool XMLTools::createFile( std::string masterTag)
{
	   // Initialize the XML4C2 system.
    try
    {
        XMLPlatformUtils::Initialize();
    }

    catch(const XMLException& toCatch)
    {
        char *pMsg = XMLString::transcode(toCatch.getMessage());
        XERCES_STD_QUALIFIER cerr << "Error during Xerces-c Initialization.\n"
             << "  Exception message:"
             << pMsg;
        XMLString::release(&pMsg);
        return false;
    }

	static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
	DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(gLS);

       if (impl != NULL)
       {
           try
           {
               m_doc = impl->createDocument(
                           0,                    // root element namespace URI.
						   XMLString::transcode(masterTag.c_str()),         // root element name
                           0);                   // document type object (DTD).

               //DOMElement* rootElem = doc->getDocumentElement();
			   }
		catch(...)
		{
			return false;
		}
		return true;
	   }
	   return false;
}