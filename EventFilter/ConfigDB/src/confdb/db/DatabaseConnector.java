package confdb.db;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.sql.SQLException;

import java.util.ArrayList;


/**
 * DatabaseConnector
 * -----------------
 * @author Philipp Schieferdecker
 * 
 * handle the low-level details of any database access, leave the DBMS
 * related details to the respective implementations of this class.
 */

abstract public class DatabaseConnector implements IDatabaseConnector
{
    //
    // member data
    //
    
    /** connection to the database */
    protected Connection connection = null;
    
    /** URL of the database */
    protected String dbURL = null;

    /** current user accessing the database */
    protected String dbUser = null;

    /** password of current user accessing the database */
    protected String dbPassword = null;
    
    //
    // construction
    //
    
    /** standard constructor */
    public DatabaseConnector(String url, String user, String password)
	throws DatabaseException
    {
	// assign the database parameters
	dbURL      = url;
	dbUser     = user;
	dbPassword = password;
	
	// open the database connection
	this.openConnection();
    }
    
    
    //
    // member functions
    //

    /** open connection to databse */
    public abstract void openConnection() throws DatabaseException;
    
    /** close the connection to the database */
    public void closeConnection() throws DatabaseException
    {
	try {
	    if (connection != null || connection.isClosed()) {
		connection.close();
		connection = null;
	    }
	}
	catch (SQLException e) {
	    String msg = "Failed to close database connection: " + e.getMessage();
	    throw new DatabaseException(msg);
	}
    }

    /** access to the connection to the database */
    public Connection getConnection()
    {
	return connection;
    }

    /** access to the URL of the database */
    public String dbURL()
    {
	return dbURL;
    }

    /** access the current username accessing the database */
    public String dbUser()
    {
	return dbUser;
    }

    /** access the password of the current user accessing the database */
    public String dbPassword()
    {
	return dbPassword;
    }
 
    /** initialize database connection parameters */
    public void setParameters(String url,String user,String password)
    {
	dbURL      = url;
	dbUser     = user;
	dbPassword = password;
    }

    /** release the resources associated with a result set */
    public ResultSet release(ResultSet rs)
    {
	if (rs != null) {
	    Statement stmt = null;
	    try {
		stmt = rs.getStatement();
	    }
	    catch (SQLException e) {
		String msg =
		    "WARNING: " +
		    "Can't retrieve SQL Statement from SQL ResultSet: " +
		    e.getMessage();
		System.out.println(msg);
	    }
	    try {
		rs.close();
	    }
	    catch (SQLException e) {
		String msg =
		    "WARNING: Can't close SQL ResultSet: "+e.getMessage();
		System.out.println(msg);
	    }
	    rs = null;
	    //if (stmt != null) {
	    //try {
	    //    stmt.close();
	    //} 
	    //catch (SQLException e) {
	    //    String msg =
	    //	"WARNING: Can't close SQL Statement: "+e.getMessage();
	    //    System.out.println(msg);
	    //}
	    //stmt = null;
	    //}
	}
	return rs;
    }

}
