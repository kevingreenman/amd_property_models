import datetime
import requests
import sys
from utils import (
    decode_and_check_mongodb_query_response,
    get_mongodb_connection_info,
    run_mongodb_query,
)

"""
After running this script, the model file can be downloaded from Dropbox as follows:

wget -O <DESIRED_FILENAME_AFTER_DOWNLOAD> <DROPBOX_SHAREABLE_LINK>

For example,

wget -O model.tar.gz https://www.dropbox.com/s/bugnhc2t6ghpucc/model.tar.gz
"""


def update_model_with_shareable_link(previous_doc, next_doc):
    http_address, mongo_port, username, password = get_mongodb_connection_info()

    current = datetime.datetime.now()
    print("Start time: " + current.strftime("%H:%M:%S.%f"))

    incoming_request = {
        "auth": {"port": mongo_port, "user": username, "password": password},
        "request": {
            "request_type": "update_model_information",
            "collection": "model_information",
            "previous_document": previous_doc,
            "new_document": next_doc,
        },
    }

    resp = requests.post(
        http_address,
        json=incoming_request,
        headers={"content-type": "application/authjson"},
    )

    current = datetime.datetime.now()
    print("End time: " + current.strftime("%H:%M:%S.%f"))

    content = decode_and_check_mongodb_query_response(resp)

    return content


def main(shareable_link, model_type="uv_vis"):
    if model_type == "uv_vis":
        mongodb_id = "629a5894199fc475e9aa3a86"
    elif model_type == "logp":
        mongodb_id = "633db3fd96f2ed8764752b44"
    elif model_type == "photodeg":
        mongodb_id = "633db3ef96f2ed8764752b43"

    content = run_mongodb_query(
        search_type="model_information",
        collection="model_information",
        search_term=["model_type", model_type],
    )

    previous_doc = content[1]["model_information_entry"][mongodb_id]
    next_doc = previous_doc.copy()
    next_doc["model_location"] = shareable_link

    content = update_model_with_shareable_link(previous_doc, next_doc)
    return


if __name__ == "__main__":
    # arguments should look like "> Share link: https://www.dropbox.com/s/u4r7lq7wpd4658v/20221020_new_uvvis_measurements_close_loop_combined_v3.tar.gz?dl=0"
    shareable_link = sys.argv[4]
    if not shareable_link.startswith("https://www.dropbox.com/s/"):
        raise ValueError("Link does not look like a Dropbox shareable link.")
    model_type = sys.argv[5]
    if model_type not in ["uv_vis", "logp", "photodeg"]:
        raise ValueError("Model type must be one of uv_vis, logp, photodeg")
    # campaign_name = sys.argv[6]
    main(shareable_link, model_type)
