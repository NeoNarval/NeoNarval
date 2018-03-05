import os
import message

from slackclient import SlackClient

authed_teams = {}

class Script(object):

    def __init__(self):
        super(Script, self).__init__()
        self.name = "NarvalSimu"
        self.emoji = ":robotface:"

        self.oauth = {"client_id": "160796522983.194115916214",
                      "client_secret": "10879f8fcb356a3942c3a96a433335e0",
                      "scope": "commands"}
        self.verification = ""
        self.client = SlackClient("")
        self.messages = {}


    def auth(self, code):

        auth_response = self.client.api_call("oauth.access",
                                            client_id=self.oauth["client_id"],
                                            client_secret=self.oauth["client_secret"],
                                            code=code
                                            )
        team_id = auth_response["team_id"]
        authed_teams[team_id] = {"bot_token":
                                auth_response["bot"]["bot_access_token"]}

        self.client = SlackClient(authed_teams[team_id]["bot_token"])