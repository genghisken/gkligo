import os, sys
import signal
from hop import Stream
from hop.auth import Auth
from hop.io import StartPosition

def handler(signum, frame):
    raise TimeoutError

class hop_reader():
    def __init__(self, username, password, topic, my_group_id, is_gcn=False):
        self.is_gcn = is_gcn
        hop_auth = Auth(username, password)
        stream   = Stream(auth=hop_auth, start_at=StartPosition.EARLIEST)
        url = f"kafka://kafka.scimma.org/{topic}"
        group_id = f"{username}-{my_group_id}"
        self.hop_stream = stream.open(url, "r", group_id=group_id).read()

    def poll(self):
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(10) 
        alert = next(self.hop_stream)
        signal.alarm(0)
        if self.is_gcn:
            return alert.fields
        else:
            return alert.content
