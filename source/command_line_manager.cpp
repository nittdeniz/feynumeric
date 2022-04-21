#include "format.hpp"
#include "messages.hpp"
#include "command_line_manager.hpp"

namespace Feynumeric
{

	Command_Line_Manager::Command_Line_Manager(int argc, char** argv){
		for( int i = 1; i < argc; ++i ){
			std::string arg{argv[i]};
//			status(FORMAT("arg: {}", arg));
			bool hasValue = false;
			for( auto j = arg.begin(); j != arg.end(); ++j ){
				if( *j == '=' ){
//					status(FORMAT("{}", std::string(arg.begin(), j)));
//					status(FORMAT("{}", std::string(j+1, arg.end())));
					arguments[std::string(arg.begin(), j)] = std::string(j + 1, arg.end());
//					status(FORMAT("{}", arguments[std::string(arg.begin(), j)]));
					hasValue = true;
					break;
				}
			}
			if( !hasValue ){
				arguments[arg] = "1";
			}
		}
	}

	void
	Command_Line_Manager::register_command(const std::string& cmd, bool is_mandatory, const std::string& description){
		commands.emplace_back(cmd, "", description, is_mandatory);
	}

	void Command_Line_Manager::register_command(std::string const& cmd, std::string const& default_value,
	                                            std::string const& description){
		commands.emplace_back(cmd, default_value, description, false);

	}

	bool Command_Line_Manager::is_enabled(const std::string& option){
		return arguments[option] == "1";
	}

	std::string Command_Line_Manager::as_string(const std::string& arg) const{
		return arguments.at(arg);
	}

	int Command_Line_Manager::as_int(const std::string& arg) const{
		return std::stoi(arguments.at(arg));
	}

	float Command_Line_Manager::as_float(const std::string& arg) const{
		return std::stof(arguments.at(arg));
	}

	double Command_Line_Manager::as_double(std::string const& arg) const{
		return std::stod(arguments.at(arg));
	}

	void Command_Line_Manager::crash_on_missing_mandatory_command() const{
		for( auto const& command : commands ){
			if( command.is_mandatory && !arguments.contains(command.cmd) ){
				for( auto const& arg : arguments ){
					status(FORMAT("Argument: {} = {}", arg.first,arg.second));
				}
				critical_error(FORMAT("Command line argument {} is mandatory but missing.", command.cmd));
			}
		}
	}
}