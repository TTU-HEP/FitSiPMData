#ifndef NTREXCEPTION_H
#define NTREXCEPTION_H

#include <string>

class NTRException
{
private:
    int line_;
    std::string funcion_;
    std::string file_;
    std::string message_;

public:
    
    NTRException(const int, const std::string, const std::string, const std::string);

    int getLineNumber() const { return line_; }
    std::string getFunctionName() const { return funcion_; }
    std::string getFileName() const { return file_; }
    std::string getMessage() const { return message_; }

    void print() const;
};

std::ostream& operator<<(std::ostream& out, const NTRException& e);

#define THROW_NTREXCEPTION( _message ) \
    throw NTRException( __LINE__, __func__, __FILE__, _message)

#endif
