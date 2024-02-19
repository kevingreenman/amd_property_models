import datetime
import json
import requests


def get_date(date_str):
    date_list = list(map(int, date_str.split("/")))
    date = datetime.date(date_list[2], date_list[0], date_list[1])
    return date


def exclude_all_but_most_recent_date(df):
    """Exclude all but the most recent date of experiments."""
    return df[df["experimentdate"] == list(df["experimentdate"])[-1]]


def get_mongodb_connection_info():
    """Get the connection information for the MongoDB database."""
    http_address = ""
    mongo_port = ""
    username = ""
    password = ""
    return http_address, mongo_port, username, password


def decode_and_check_mongodb_query_response(resp):
    """Decode the response from the MongoDB query and check for errors.

    Args:
        resp: Response from the MongoDB query.

    Returns:
        content: Decoded response from the MongoDB query.
    """
    print("Respose status code: " + str(resp.status_code))
    if resp.status_code != 200:
        raise RuntimeError("Something went wrong with your request")
    else:
        try:
            content = json.loads(resp.content.decode("utf_8"))
        except json.JSONDecodeError as e:
            print(e)
            raise RuntimeError("Something went wrong with decoding the content JSON...")
        try:
            print("Content print statements: ", content[1]["print_statements"])
        except KeyError as e:
            print(e)
            print("Content data errors: ", content[1]["data_errors"])
            raise RuntimeError(
                "Something went wrong with printing the print statement..."
            )
        return content


def run_mongodb_query(
    search_type: str,
    collection: str,
    search_term: list,
):
    """Run a MongoDB query.

    Args:
        search_type (str): Type of search to run. Options are "special_search" and "model_information".
        collection (str): Collection to search. Options are "data" and "model_information".
        search_term (list): List of search terms. Options include but are not limited to ["exp_abs_spec"], ['exp_abs_stick'], ['model_type', 'uv_vis'], and ["exp_pl_spec"]

    Returns:
        dict: Decoded response from the MongoDB query.
    """
    http_address, mongo_port, username, password = get_mongodb_connection_info()

    current = datetime.datetime.now()
    print("Start time: " + current.strftime("%H:%M:%S.%f"))

    incoming_request = {
        "auth": {"port": mongo_port, "user": username, "password": password},
        "request": {
            "request_type": "query_data",
            "search_type": search_type,
            "collection": collection,
            "document_property": "data_entry.type_tag",
            "search_term": search_term,
        },
    }

    resp = requests.post(
        http_address,
        json=incoming_request,
        headers={"content-type": "application/authjson"},
    )

    content = decode_and_check_mongodb_query_response(resp)

    current = datetime.datetime.now()
    print("End time: " + current.strftime("%H:%M:%S.%f"))

    return content
